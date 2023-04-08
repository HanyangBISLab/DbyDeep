import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import pandas as pd
import tensorflow as tf
import os

from sklearn.model_selection import train_test_split
import numpy as np
import argparse


def get_npy_DbyDeep(df):
    label_enc = {v:k for k, v in enumerate('ZARNDCQEGHILKMFPSTWYV')}  # Z : 0
    pep_data = [[label_enc[aa] for aa in seq] + [0]*(40-len(seq))  # zero padding
               for seq in df.peptide.values]
    nterm_data = [[label_enc[aa] for aa in seq]
               for seq in df.nterm.values]
    cterm_data = [[label_enc[aa] for aa in seq]
               for seq in df.cterm.values]
    miss1_data = [[label_enc[aa] for aa in seq]
               for seq in df.miss1.values]
    miss2_data = [[label_enc[aa] for aa in seq]
               for seq in df.miss2.values]
    return np.array(pep_data), np.array(nterm_data), np.array(cterm_data), np.array(miss1_data), np.array(miss2_data), np.array(df.label.values)

def main(opt):
    # data load
    print("data loading...")
    data = pd.read_csv(opt.data_path)
    
    # model load
    print("model loading...")
    model = tf.keras.models.load_model(opt.model_path)

    # prediction
    if opt.retrain_flag:
        print('retraining...')
        df_test=data
        res=df_test.label.value_counts()  # downsampling
        lab_cnts={k:v for k,v in zip(res.index, res.values)}
        num=min(lab_cnts.values())
        _=df_test.loc[df_test.label==True].sample(num, random_state=2023)
        __=df_test.loc[df_test.label==False].sample(num, random_state=2023)
        df_test=pd.concat([_,__],axis=0).reset_index(drop=True)
        df_test = df_test[df_test.peptide.apply(lambda x: 'B' not in x)].reset_index(drop=True)
        df_test = df_test[df_test.nterm.apply(lambda x: 'B' not in x)].reset_index(drop=True)
        df_test = df_test[df_test.cterm.apply(lambda x: 'B' not in x)].reset_index(drop=True)
        df_test = df_test[df_test.miss1.apply(lambda x: 'B' not in x)].reset_index(drop=True)
        df_test = df_test[df_test.miss2.apply(lambda x: 'B' not in x)].reset_index(drop=True)
        df_test = df_test[df_test.peptide.apply(lambda x: 'Z' not in x)].reset_index(drop=True)
        data=df_test
        data_train, data_val = train_test_split(data, test_size=0.1, random_state=2023)
        pep_train, n_train, c_train, m1_train, m2_train, label_train = get_npy_DbyDeep(data_train)
        pep_val, n_val, c_val, m1_val, m2_val, label_val = get_npy_DbyDeep(data_val)
        model.compile(loss=tf.keras.losses.BinaryCrossentropy(from_logits=True),
                optimizer=tf.keras.optimizers.Adam(1e-5),  # low learning rate
                metrics=['accuracy'])
        es = tf.keras.callbacks.EarlyStopping(monitor='val_loss',
                                            mode='min', 
                                            verbose=1,
                                            patience=50)
        cp = tf.keras.callbacks.ModelCheckpoint(
            f'{opt.save_path}DbyDeep_retrained.h5',
            monitor='val_loss', 
            verbose=0, 
            save_best_only=True, 
            save_weights_only=False, 
            mode='auto', 
            period=1)
        model.fit([pep_train, n_train, c_train, m1_train, m2_train], label_train,
                    epochs=20,
                    batch_size=256,
                    validation_data=([pep_val, n_val, c_val, m1_val, m2_val], label_val),
                    callbacks=[es, cp],
                    )
    else:
        print("prediction...")
        pep_test, n_test, c_test, m1_test, m2_test, label_test = get_npy_DbyDeep(data)
        probs = model.predict([pep_test, n_test, c_test, m1_test, m2_test])
        y_pred = [1 if i>=0.5 else 0 for i in probs]
        data['Prob'] = probs
        data['Detectability'] = y_pred
        data.to_csv(opt.save_path+opt.job_name+'.csv', index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--retrain-flag', type=bool, default=False, help='retrain flag')
    parser.add_argument('--data-path', type=str, default='./data/test.csv', help='data-path')
    parser.add_argument('--model-path', type=str, default='./data/DbyDeep.h5', help='model-path')
    parser.add_argument('--save-path', type=str, default='./data/', help='save-path')
    parser.add_argument('--job-name', type=str, default='test_result', help='save-path')
    opt = parser.parse_args()
    if opt.save_path[-1]!='/':
        opt.save_path+='/'
    print(opt)
    main(opt)
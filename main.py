import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import pandas as pd
import tensorflow as tf
import numpy as np


def get_npy_DbyDeep(df):
    label_enc = {v:k for k, v in enumerate('ZARNDCQEGHILKMFPSTWYV')}  # Z : 0
    pep_data = [[label_enc[aa] for aa in seq] + [0]*(81-len(seq))  # zero padding
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


def main():
    # gpu setting
    # gpus = tf.config.experimental.list_physical_devices('GPU')
    # if gpus:
    #     try:
    #         tf.config.experimental.set_virtual_device_configuration(gpus[0],
    #         [tf.config.experimental.VirtualDeviceConfiguration(memory_limit=1024*4)])
    #     except RuntimeError as e:
    #         print(e)
    
    
    # data load
    print("data loading...")
    data = pd.read_csv('data/test.csv')
    pep_test, n_test, c_test, m1_test, m2_test, label_test = get_npy_DbyDeep(data)
    
    # model load
    print("model loading...")
    model = tf.keras.models.load_model("data/DbyDeep.h5")

    # prediction
    print("prediction...")
    probs = model.predict([pep_test, n_test, c_test, m1_test, m2_test])
    y_pred = [1 if i>=0.5 else 0 for i in probs]
    data['Prob'] = probs
    data['Detectability'] = y_pred
    data.to_csv('data/test_result.csv', index=False)


if __name__ == "__main__":
    main()
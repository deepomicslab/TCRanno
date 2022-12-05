from numpy import *
import os,sys
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'
import tensorflow as tf
from tensorflow.keras.models import Sequential, Model
from tensorflow.keras.layers import Layer, Dense, Activation, Conv1D, Flatten, Dropout, Input, BatchNormalization, Reshape
from tensorflow.keras import backend as K
from tensorflow import keras

amino_acids=['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V','-']
def Onehot_encoding(seqs, aa=amino_acids, max_length=30):
    n=len(seqs)
    matrix=zeros((n,max_length,21),dtype=int8)
    for x in range(n):
        seq=seqs[x]+(max_length-len(seqs[x]))*'-'
        for i in range(max_length):
            matrix[x][i][aa.index(seq[i].upper())] = 1
    return matrix

def get_norm_latent(inseqs,encoder):
    X = Onehot_encoding(inseqs)
    latent = encoder.predict(X)
    norm_latent = array([j/sqrt(sum(j**2)) for j in latent])
    return norm_latent

class dSampling(Layer): # dummy sampling layer, epsilon set to 0 so that prediction results are repeatable
    def call(self, inputs):
        z_mean, z_log_var = inputs
        #batch = tf.shape(z_mean)[0]
        #dim = tf.shape(z_mean)[1]
        epsilon = 0 #epsilon = tf.keras.backend.random_normal(shape=(batch, dim),mean=0,stddev=1e-12)
        return z_mean + tf.exp(0.5 * z_log_var) * epsilon

def load_encoder(model_path=None):
    if model_path==None:
        d = os.path.dirname(sys.modules['tcranno'].__file__)
        resouce = 'pretrained/pretrained_encoder.h5'
        model_path = os.path.join(d,resouce)
    encoder = keras.models.load_model(model_path, compile=False, custom_objects={'Sampling1':dSampling})
    return encoder
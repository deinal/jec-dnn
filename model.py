import tensorflow as tf
from tensorflow.keras import Input, Model
from tensorflow.keras.layers import Activation, Dense, TimeDistributed, Layer
from tensorflow.keras.utils import plot_model


class Sum(Layer):
    def call(self, inputs: tf.RaggedTensor) -> tf.Tensor:
        return tf.math.reduce_sum(inputs, axis=1)


def get_model(num_features):
    inputs = Input(shape=(None, num_features), ragged=True, name='inputs')

    inputs_slice = Input(shape=(inputs.shape[-1],), name='inputs_slice')

    x = Dense(100, kernel_initializer='he_uniform', use_bias=True, name='deepset_dense')(inputs_slice)

    outputs_slice = Activation('relu', name='deepset_activation')(x)

    submodel_slice = Model(inputs=inputs_slice, outputs=outputs_slice, name='deepset_submodel')

    outputs = TimeDistributed(submodel_slice, name='deepset_distributed')(inputs)

    outputs = Sum(name='deepset_sum')(outputs)

    x = Dense(100, kernel_initializer='he_uniform', use_bias=True, name='head_dense')(outputs)

    x = Activation('relu', name='head_activation')(x)

    outputs = Dense(1, name='head_dense_output')(x)

    model = Model(inputs=inputs, outputs=outputs, name='full')

    _summarize_model(model)

    return model


def _summarize_model(model):
    for layer in model.layers:
        if isinstance(layer, TimeDistributed):
            layer.layer.summary()

    model.summary()

    plot_model(model, 'plots/model.png', dpi=100, show_shapes=True, expand_nested=True)

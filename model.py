import tensorflow as tf
from tensorflow.keras import Input, Model
from tensorflow.keras.layers import Activation, Dense, TimeDistributed, BatchNormalization, Dropout,  Concatenate, Layer


class Sum(Layer):
    def __init__(self, axis=1, **kwargs):
        super().__init__(**kwargs)
        self.axis = axis

    def call(self, inputs):
        return tf.math.reduce_sum(inputs, axis=self.axis)


def get_model(num_features, num_globals, config):
    constituents = Input(shape=(None, num_features), ragged=True, name='constituents')

    constituents_slice = Input(shape=(constituents.shape[-1],), name='constituents_slice')

    deepset_outputs_slice = _mlp(constituents_slice, config, 'deepset')

    deepset_model_slice = Model(inputs=constituents_slice, outputs=deepset_outputs_slice, name='deepset_model_slice')

    deepset_outputs = TimeDistributed(deepset_model_slice, name='deepset_distributed')(constituents)

    constituents_head = Sum(name='constituents_head')(deepset_outputs)

    globals = Input(shape=(num_globals,), name='globals')

    inputs_head = Concatenate(name='head')([constituents_head] + [globals])

    x = _mlp(inputs_head, config, 'head')

    outputs = Dense(1, name='head_dense_output')(x)

    model = Model(inputs=[constituents, globals], outputs=outputs, name='dnn')

    model.summary()

    for layer in model.layers:
        if isinstance(layer, TimeDistributed):
            layer.layer.summary()

    return model


def _mlp(x, config, name):
    for i, units in enumerate(config['layers'], start=1):
        x = Dense(units, kernel_initializer=config['initializer'], name=f'{name}_dense_{i}')(x)
        x = BatchNormalization(name=f'{name}_batch_normalization_{i}')(x)
        x = Activation(config['activation'], name=f'{name}_activation_{i}')(x)
        x = Dropout(config['dropout'], name=f'{name}_dropout_{i}')(x)
    return x

import tensorflow as tf
from tensorflow.keras import Input, Model
from tensorflow.keras.layers import Activation, Dense, TimeDistributed, Layer
from tensorflow.keras.utils import plot_model


class Sum(Layer):
    def call(self, inputs: tf.RaggedTensor) -> tf.Tensor:
        return tf.math.reduce_sum(inputs, axis=1)


def get_model(outdir, num_features, config):
    inputs = Input(shape=(None, num_features), ragged=True, name='inputs')

    inputs_slice = Input(shape=(inputs.shape[-1],), name='inputs_slice')

    outputs_slice = _mlp(inputs_slice, config['layers'], config['activation'], config['initializer'], 'deepset')

    deepset_model_slice = Model(inputs=inputs_slice, outputs=outputs_slice, name='deepset_model_slice')

    deepset_outputs = TimeDistributed(deepset_model_slice, name='deepset_distributed')(inputs)

    head_inputs_deepset = Sum(name='deepset_sum')(deepset_outputs)

    x = _mlp(head_inputs_deepset, config['layers'], config['activation'], config['initializer'], 'head')

    outputs = Dense(1, name='head_dense_output')(x)

    model = Model(inputs=inputs, outputs=outputs, name='full')

    _summarize_model(outdir, model)

    return model


def _mlp(x, layers, activation, initializer, name):
    for i, units in enumerate(layers):
        x = Dense(units, kernel_initializer=initializer, use_bias=True, name=f'{name}_dense_{i}')(x)
        x = Activation(activation, name=f'{name}_activation_{i}')(x)
    return x


def _summarize_model(outdir, model):
    for layer in model.layers:
        if isinstance(layer, TimeDistributed):
            layer.layer.summary()

    model.summary()

    plot_model(model, f'./results/{outdir}/model.pdf', dpi=100, show_shapes=True, expand_nested=True)

import tensorflow as tf
from tensorflow.keras import Input, Model
from tensorflow.keras.layers import Activation, Dense, TimeDistributed, BatchNormalization, Dropout,  Concatenate, Layer
from tensorflow.keras.utils import plot_model


class Sum(Layer):
    def call(self, inputs: tf.RaggedTensor) -> tf.Tensor:
        return tf.math.reduce_sum(inputs, axis=1)


def get_model(outdir, num_features, num_globals, config):
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

    _summarize_model(outdir, model)

    return model


def _mlp(x, config, name):
    for i, units in enumerate(config['layers'], start=1):
        x = Dense(units, kernel_initializer=config['initializer'], use_bias=True, name=f'{name}_dense_{i}')(x)
        x = BatchNormalization(name=f'{name}_batch_normalization_{i}')(x)
        x = Activation(config['activation'], name=f'{name}_activation_{i}')(x)
        x = Dropout(config['dropout'], name=f'{name}_dropout_{i}')(x)
    return x


def _summarize_model(outdir, model):
    for layer in model.layers:
        if isinstance(layer, TimeDistributed):
            layer.layer.summary()
            plot_model(layer.layer, f'./results/{outdir}/{layer.layer.name}.pdf', dpi=100, show_shapes=True)

    model.summary()

    plot_model(model, f'./results/{outdir}/model.pdf', dpi=100, show_shapes=True)
    plot_model(model, f'./results/{outdir}/full_model.pdf', dpi=100, show_shapes=True, expand_nested=True)

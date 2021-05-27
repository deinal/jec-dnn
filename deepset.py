from tensorflow.keras import Input, Model
from tensorflow.keras.layers import Activation, Dense, TimeDistributed, BatchNormalization, Dropout, Concatenate, Add
from layers import Sum


def get_deepset(num_constituents, num_globals, config):
    constituents = Input(shape=(None, num_constituents), ragged=True, name='constituents')

    constituents_slice = Input(shape=(constituents.shape[-1],), name='constituents_slice')

    if config['type'] == 'mlp':
        deepset_outputs_slice = _mlp(constituents_slice, config, name='deepset')
    if config['type'] == 'resnet':
        deepset_outputs_slice = _resnet(constituents_slice, config, name='deepset')

    deepset_model_slice = Model(inputs=constituents_slice, outputs=deepset_outputs_slice, name='deepset_model_slice')

    deepset_outputs = TimeDistributed(deepset_model_slice, name='deepset_distributed')(constituents)

    constituents_head = Sum(axis=1, name='constituents_head')(deepset_outputs)

    globals = Input(shape=(num_globals,), name='globals')

    inputs_head = Concatenate(name='head')([constituents_head, globals])

    if config['type'] == 'mlp':
        x = _mlp(inputs_head, config, name='head')
    if config['type'] == 'resnet':
        x = _resnet(inputs_head, config, name='head')

    outputs = Dense(1, name='head_dense_output')(x)

    model = Model(inputs=[constituents, globals], outputs=outputs, name='dnn')

    model.summary()

    for layer in model.layers:
        if isinstance(layer, TimeDistributed):
            layer.layer.summary()

    return model


def _mlp(x, config, name):
    for idx, units in enumerate(config['units'], start=1):
        x = Dense(units, kernel_initializer=config['initializer'], name=f'{name}_dense_{idx}')(x)
        x = BatchNormalization(name=f'{name}_batch_normalization_{idx}')(x)
        x = Activation(config['activation'], name=f'{name}_activation_{idx}')(x)
        x = Dropout(config['dropout'], name=f'{name}_dropout_{idx}')(x)
    return x


def _resnet(x, config, name):
    units = config['units']
    for idx in range(0, len(units) - 1, 2):
        n1 = units[idx]
        n2 = units[idx + 1]
        layer_idx = idx // 2 + 1

        x_main = Dense(n1, kernel_initializer=config['initializer'], name=f'{name}_dense_{layer_idx}_1')(x)
        x_main = Activation(config['activation'], name=f'{name}_activation_{layer_idx}_1')(x_main)
        x_main = Dense(n2, kernel_initializer=config['initializer'], name=f'{name}_dense_{layer_idx}_2')(x_main)

        # Include a projection to match the dimensions
        sc = Dense(n2, kernel_initializer=config['initializer'], use_bias=False, name=f'{name}_projection_{layer_idx}')(x)

        x = Add(name=f'add_{layer_idx}')([x_main, sc])
        x = Activation(config['activation'], name=f'{name}_activation_{layer_idx}_2')(x)

    if len(units) % 2 == 1:
        x = Dense(units[-1], kernel_initializer=config['initializer'], name=f'{name}_dense_last')(x)
        x = Activation(config['activation'], name=f'{name}_activation_last')(x)

    return x

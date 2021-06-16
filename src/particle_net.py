import tensorflow as tf
from tensorflow.keras import Input, Model
from tensorflow.keras.layers import Activation, Add, BatchNormalization, Conv2D, Dense, Dropout, Layer, Multiply, Concatenate
from src.layers import Mean, Max, Expand, Squeeze


def get_particle_net(num_constituents, num_points, num_globals, config):
    """
    ParticleNet: Jet Tagging via Particle Clouds
    arxiv.org/abs/1902.08570
    
    Parameters
    ----------
    input_shapes : dict
        The shapes of each input (`points`, `features`, `mask`).
    """

    points = Input(name='points', shape=(num_points, 2))
    features = Input(name='features', shape=(num_points, num_constituents))
    mask = Input(name='mask', shape=(num_points, 1))
    coord_shift = Input(name='coord_shift', shape=(num_points, 1))
    globals = Input(name='globals', shape=(num_globals,))

    outputs = _particle_net_base(points, features, mask, coord_shift, globals, config)

    model = Model(inputs=[points, features, mask, coord_shift, globals], outputs=outputs)

    model.summary()

    return model


def _particle_net_base(points, features, mask, coord_shift, globals, config):
    """
    points : (N, P, C_coord)
    features:  (N, P, C_features), optional
    mask: (N, P, 1), optional
    """

    # fts = tf.squeeze(BatchNormalization(name='fts_bn')(tf.expand_dims(features, axis=2)), axis=2)
    fts = features
    for layer_idx, channels in enumerate(config['channels'], start=1):
        pts = Add(name=f'add_{layer_idx}')([coord_shift, points]) if layer_idx == 1 else Add(name=f'add_{layer_idx}')([coord_shift, fts])
        fts = _edge_conv(
            pts, fts, config['num_points'], channels, config, name=f'edge_conv_{layer_idx}'
        )

    fts = Multiply()([fts, mask])

    pool = Mean(axis=1)(fts) # (N, C)

    x = Concatenate(name='head')([pool, globals])

    # x = pool
    for layer_idx, units in enumerate(config['units']):
        x = Dense(units)(x)
        x = Activation(config['activation'])(x)
        if config['dropout']:
            x = Dropout(config['dropout'])(x)
    out = Dense(1, name='out')(x)
    return out # (N, num_classes)


def _edge_conv(points, features, num_points, channels, config, name):
    """EdgeConv
    Args:
        K: int, number of neighbors
        in_channels: # of input channels
        channels: tuple of output channels
        pooling: pooling method ('max' or 'average')
    Inputs:
        points: (N, P, C_p)
        features: (N, P, C_0)
    Returns:
        transformed points: (N, P, C_out), C_out = channels[-1]
    """

    fts = features
    knn_fts = KNearestNeighbors(num_points, config['K'], name=f'{name}_knn')([points, fts])

    x = knn_fts
    for idx, channel in enumerate(channels, start=1):
        x = Conv2D(
            channel, kernel_size=(1, 1), strides=1, data_format='channels_last',
            use_bias=False if config['batch_norm'] else True, kernel_initializer=config['initializer'], name=f'{name}_conv_{idx}'
        )(x)
        if config['batch_norm']:
            x = BatchNormalization(name=f'{name}_batchnorm_{idx}')(x)
        if config['activation']:
            x = Activation(config['activation'], name=f'{name}_activation_{idx}')(x)

    if config['pooling'] == 'max':
        fts = Max(axis=2, name=f'{name}_max')(x) # (N, P, C')
    else:
        fts = Mean(axis=2, name=f'{name}_mean')(x) # (N, P, C')

    if config['shortcut']:
        sc = Expand(axis=2, name=f'{name}_shortcut_expand')(features)
        sc = Conv2D(
            channels[-1], kernel_size=(1, 1), strides=1, data_format='channels_last',
            use_bias=False if config['batch_norm'] else True, kernel_initializer=config['initializer'], name=f'{name}_shortcut_conv'
        )(sc)
        if config['batch_norm']:
            sc = BatchNormalization(name=f'{name}_shortcut_batchnorm')(sc)
        sc = Squeeze(axis=2, name=f'{name}_shortcut_squeeze')(sc)

        x = Add(name=f'{name}_add')([sc, fts])
    else:
        x = fts

    return Activation(config['activation'], name=f'{name}_activation')(x) # (N, P, C')


class KNearestNeighbors(Layer):
    def __init__(self, num_points, k, **kwargs):
        super().__init__(**kwargs)
        self.num_points = num_points
        self.k = k

    def call(self, inputs):
        points, features = inputs
        # distance
        D = _batch_distance_matrix_general(points, points) # (N, P, P)
        _, top_k_indices = tf.math.top_k(-D, k=self.k + 1) # (N, P, K+1)
        top_k_indices = top_k_indices[:, :, 1:] # (N, P, K)

        queries_shape = tf.shape(features)
        batch_size = queries_shape[0]
        batch_indices = tf.tile(tf.reshape(tf.range(batch_size), (-1, 1, 1, 1)), (1, self.num_points, self.k, 1))
        indices = tf.concat([batch_indices, tf.expand_dims(top_k_indices, axis=3)], axis=3) # (N, P, K, 2)
        
        knn_fts =  tf.gather_nd(features, indices) # (N, P, K, C)
        knn_fts_center = tf.tile(tf.expand_dims(features, axis=2), (1, 1, self.k, 1)) # (N, P, K, C)

        return tf.concat([knn_fts_center, tf.subtract(knn_fts, knn_fts_center)], axis=-1) # (N, P, K, 2*C)


# A shape is (N, P_A, C), B shape is (N, P_B, C)
# D shape is (N, P_A, P_B)
def _batch_distance_matrix_general(A, B):
    r_A = tf.math.reduce_sum(A * A, axis=2, keepdims=True)
    r_B = tf.math.reduce_sum(B * B, axis=2, keepdims=True)
    m = tf.linalg.matmul(A, tf.transpose(B, perm=(0, 2, 1)))
    D = r_A - 2 * m + tf.transpose(r_B, perm=(0, 2, 1))
    return D

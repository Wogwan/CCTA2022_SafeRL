import tensorflow as tf
import numpy as np
from utils import *
from sklearn.linear_model import LinearRegression

# Barrier Function Compensator
class BARRIER():
    def __init__(self, input_size, action_size, learning_rate, mlp_mode):
        self.input_size = input_size
        self.action_size = action_size
        self.learning_rate = learning_rate
        self.build_model()
        if mlp_mode == 'original':
            # self.train = self.train_revised
            self.train = self.train_original

    def create_mlp_model(self):
        # Input()
        inputs = tf.keras.Input(shape=(self.input_size,), name='Obs')
        out = tf.keras.layers.Dense(units=32, use_bias=True, name='h1', activation='relu',
                                    kernel_initializer=tf.keras.initializers.TruncatedNormal(stddev=0.01))(inputs)
        out = tf.keras.layers.Dense(units=8, use_bias=True, name='h2', activation='relu',
                                    kernel_initializer=tf.keras.initializers.TruncatedNormal(stddev=0.01))(out)
        out = tf.keras.layers.Dense(units=self.action_size, use_bias=True, name='h3',
                                    kernel_initializer=tf.keras.initializers.TruncatedNormal(stddev=0.01))(out)
        mlp_model = tf.keras.models.Model(inputs=inputs, outputs=out)
        return mlp_model

    # Input will be state : [batch_size, observation size]
    # Ouput will be control action
    def build_model(self):
        print('Initializing Barrier Compensation network')
        self.mlp_model = self.create_mlp_model()
        self.mlp_model.summary()
        self.mlp_opt = tf.keras.optimizers.Adam(self.learning_rate)

        # with tf.variable_scope('Compensator'):
        #     #Input will be observation
        #     self.x = tf.placeholder(tf.float32, [None, self.input_size], name='Obs')
        #     #Target will be control action
        #     self.target = tf.placeholder(tf.float32, [None, self.action_size], name='Target_comp')

        #     #Model is MLP composed of 2 hidden layers with 50, 40 relu units
        #     h1 = LINEAR(self.x, 30, name='h1')
        #     h1_n1 = tf.nn.relu(h1)
        #     h2 = LINEAR(h1_n1, 20, name='h2')
        #     h2_n1 = tf.nn.relu(h2)
        #     self.value = LINEAR(h2_n1, self.action_size, name='h3')

        # tr_vrbs = self.mlp_model.trainable_variables
        # for i in tr_vrbs:
        #     print(i.op.name)

        # Compute the loss and gradient of loss w.r.t. neural network weights
        # target: historical u_BAR + current u_bar组成当前k的CBF补偿
        # value: 输入当前k的起始状态s_t, 通过MLP去模拟当前K iter的补偿，然后计算loss.

        # To adjust weights
        # self.get_value = GetValue(self.sess, tr_vrbs, name='Compensator')
        # self.set_value = SetValue(self.sess, tr_vrbs, name='Compensator')

        # self.sess.run(tf.global_variables_initializer())

    def get_training_rollouts(self, paths):
        # Get observations and actions
        self.observation = np.squeeze(np.concatenate([path["Observation"] for path in paths]))
        self.action_bar = np.concatenate([path["Action_bar"] for path in paths])
        self.action_BAR = np.concatenate([path["Action_BAR"] for path in paths])

        # Reshape observations & actions to use for training
        batch_s = self.observation.shape[0]
        self.action_bar = np.resize(self.action_bar, [batch_s, self.action_size])
        self.action_BAR = np.resize(self.action_BAR, [batch_s, self.action_size])

    # Given current observation, get the neural network output (representing barrier compensator)
    def get_action(self, obs):
        # (3, ) => (1, 3)
        observation = np.expand_dims(np.squeeze(obs), 0)
        u_bar = self.mlp_model.predict(observation)

        # feed_dict = {self.x:observation}
        # u_bar = self.sess.run(self.value, feed_dict)
        return u_bar

    def train_revised(self):
        # print('Training barrier function compensator')
        action_comp = self.action_bar + self.action_BAR

        with tf.GradientTape() as tape:
            value_pred = self.mlp_model(self.observation, training=True)
            self.loss = tf.math.reduce_mean(tf.math.pow(action_comp - value_pred, 2))
        grads = tape.gradient(self.loss, self.mlp_model.trainable_variables)
        self.mlp_opt.apply_gradients(zip(grads, self.mlp_model.trainable_variables))
        return self.loss
        # gradient_objective = grads # tf.concat(values=[tf.reshape(g, [np.prod(v.get_shape().as_list()),]) for (g,v) in zip(grads, self.mlp_model.trainable_variables)], axis=0)

        """
        for i in range(10):
            #Get the parameter values for gradient, etc...
            parameter_prev = self.mlp_model.get_weights()
            # parameter_prev = self.get_value()
            action_comp = self.action_bar + self.action_BAR

            with tf.GradientTape() as tape:
                value_pred = self.mlp_model(self.observation, training=True)
                self.loss = tf.math.reduce_mean(tf.math.pow(action_comp - value_pred, 2))
            grads = tape.gradient(self.loss, self.mlp_model.trainable_variables)
            gradient_objective = grads # tf.concat(values=[tf.reshape(g, [np.prod(v.get_shape().as_list()),]) for (g,v) in zip(grads, self.mlp_model.trainable_variables)], axis=0)


            # gradient_objective = FLAT_GRAD(self.loss, self.mlp_model.trainable_variables) 

            #Function which takes 'y' input and returns Hy
            def get_hessian_vector_product(y):

                # self.HVP = HESSIAN_VECTOR_PRODUCT(self.loss, self.mlp_model.trainable_variables, y)
                self.HVP = HESSIAN_VECTOR_PRODUCT(self, y, action_comp)
                return self.HVP

            #Get loss under current parameter
            def loss_func(parameter):
                # self.set_value(parameter)
                self.mlp_model.set_weights(parameter)
                value_pred = self.mlp_model(self.observation, training=True)
                self.loss = tf.math.reduce_mean(tf.math.pow(action_comp - value_pred, 2))
                return self.loss
            '''
            #Move theta in direction that minimizes loss (improve barrier function parameterization)
            step_direction = CONJUGATE_GRADIENT(get_hessian_vector_product, -gradient_objective)

            #Determine step to satisfy contraint and minimize loss
            constraint_approx = 0.5*step_direction.dot(get_hessian_vector_product(step_direction))
            maximal_step_length = np.sqrt(self.args.bar_constraint_max /constraint_approx)
            full_step = maximal_step_length*step_direction
            '''
            bar_constraint_max = 0.02

            full_step = []
            for sig_element in gradient_objective:
                scalar_element = -bar_constraint_max * sig_element.numpy()
                full_step.append(scalar_element)
            # full_step = -bar_constraint_max * np.array(gradient_objective)
            #Set parameter to decrease loss - use line search to check
            new_parameter = LINE_SEARCH(loss_func, parameter_prev, full_step, name='Barrier loss')
            # self.set_value(new_parameter, update_info=0)
            self.mlp_model.set_weights(new_parameter)

            if (np.array_equal(new_parameter, parameter_prev)):
                print("Break")
                #continue
                return loss_func(new_parameter)
        return loss_func(new_parameter)
        """

    def train_original(self):
        # print('Training barrier function compensator')
        action_comp = self.action_bar + self.action_BAR

        for i in range(10):
            # Get the parameter values for gradient, etc...
            parameter_prev = self.mlp_model.get_weights()
            # parameter_prev = self.get_value()
            action_comp = self.action_bar + self.action_BAR

            with tf.GradientTape() as tape:
                value_pred = self.mlp_model(self.observation, training=True)
                self.loss = tf.math.reduce_mean(tf.math.pow(action_comp - value_pred, 2))
            grads = tape.gradient(self.loss, self.mlp_model.trainable_variables)
            gradient_objective = grads  # tf.concat(values=[tf.reshape(g, [np.prod(v.get_shape().as_list()),]) for (g,v) in zip(grads, self.mlp_model.trainable_variables)], axis=0)

            # gradient_objective = FLAT_GRAD(self.loss, self.mlp_model.trainable_variables)

            # Function which takes 'y' input and returns Hy
            def get_hessian_vector_product(y):

                # self.HVP = HESSIAN_VECTOR_PRODUCT(self.loss, self.mlp_model.trainable_variables, y)
                self.HVP = HESSIAN_VECTOR_PRODUCT(self, y, action_comp)
                return self.HVP

            # Get loss under current parameter
            def loss_func(parameter):
                # self.set_value(parameter)
                self.mlp_model.set_weights(parameter)
                value_pred = self.mlp_model(self.observation, training=True)
                self.loss = tf.math.reduce_mean(tf.math.pow(action_comp - value_pred, 2))
                return self.loss

            '''
            #Move theta in direction that minimizes loss (improve barrier function parameterization)
            step_direction = CONJUGATE_GRADIENT(get_hessian_vector_product, -gradient_objective)

            #Determine step to satisfy contraint and minimize loss
            constraint_approx = 0.5*step_direction.dot(get_hessian_vector_product(step_direction))
            maximal_step_length = np.sqrt(self.args.bar_constraint_max /constraint_approx)
            full_step = maximal_step_length*step_direction
            '''
            bar_constraint_max = 0.02

            full_step = []
            for sig_element in gradient_objective:
                scalar_element = -bar_constraint_max * sig_element.numpy()
                full_step.append(scalar_element)
            # full_step = -bar_constraint_max * np.array(gradient_objective)
            # Set parameter to decrease loss - use line search to check
            new_parameter = LINE_SEARCH(loss_func, parameter_prev, full_step, name='Barrier loss')
            # self.set_value(new_parameter, update_info=0)
            self.mlp_model.set_weights(new_parameter)

            if (np.array_equal(new_parameter, parameter_prev)):
                print("Break")
                # continue
                return loss_func(new_parameter)
        return loss_func(new_parameter)

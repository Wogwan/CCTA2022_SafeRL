"""
Last update: 2022-06-17
Name: ddpg_qp.py
Description: Implementation of MLP model for learning the previous compensation from the controller.
"""

import tensorflow as tf
import numpy as np
from utils import HESSIAN_VECTOR_PRODUCT, LINE_SEARCH

# Barrier Function Compensator
class BARRIER():
    def __init__(self, input_size, action_size, learning_rate, mlp_mode):
        self.input_size = input_size
        self.action_size = action_size
        self.learning_rate = learning_rate
        self.build_model()
        if mlp_mode == 'linesearch':
            self.train = self.train_linesearch
        elif mlp_mode == 'adam':
            self.train = self.train_optimizer
        else:
            raise Exception('Cannot find the mlp mode [{}]'.format(mlp_mode))

    def create_mlp_model(self):
        inputs = tf.keras.Input(shape=(self.input_size,), name='Obs')
        out = tf.keras.layers.Dense(units=32, use_bias=True, name='h1', activation='relu', 
                                    kernel_initializer=tf.keras.initializers.TruncatedNormal(stddev=0.01))(inputs)
        out = tf.keras.layers.Dense(units=8, use_bias=True, name='h2', activation='relu', 
                                    kernel_initializer=tf.keras.initializers.TruncatedNormal(stddev=0.01))(out)
        out = tf.keras.layers.Dense(units=self.action_size, use_bias=True, name='h3', 
                                    kernel_initializer=tf.keras.initializers.TruncatedNormal(stddev=0.01))(out)
        mlp_model = tf.keras.models.Model(inputs=inputs, outputs=out)
        return mlp_model

    def build_model(self):
        """
        Input: [batch_size, observation size]
        Ouput: control action
        """
        print('Initializing Barrier Compensation network')
        self.mlp_model = self.create_mlp_model()
        self.mlp_model.summary()
        self.mlp_opt = tf.keras.optimizers.Adam(self.learning_rate)

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
        return u_bar

    def train_optimizer(self):
        # The MLP model is trained by the adam optimizer.
        action_comp = self.action_bar + self.action_BAR
        with tf.GradientTape() as tape:
            value_pred = self.mlp_model(self.observation, training=True)
            self.loss = tf.math.reduce_mean(tf.math.pow(action_comp - value_pred, 2))
        grads = tape.gradient(self.loss, self.mlp_model.trainable_variables)
        self.mlp_opt.apply_gradients(zip(grads, self.mlp_model.trainable_variables))
        return self.loss

    def train_linesearch(self):
        # The MLP model is trained by line search method.
        action_comp = self.action_bar + self.action_BAR            
        for i in range(10):
            parameter_prev = self.mlp_model.get_weights()
            action_comp = self.action_bar + self.action_BAR

            with tf.GradientTape() as tape:
                value_pred = self.mlp_model(self.observation, training=True)
                self.loss = tf.math.reduce_mean(tf.math.pow(action_comp - value_pred, 2)) # MSE
            grads = tape.gradient(self.loss, self.mlp_model.trainable_variables)
            gradient_objective = grads
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
            grad_step = maximal_step_length*step_direction
            '''
            bar_constraint_max = 0.02
            grad_step = []
            for sig_element in gradient_objective:
                scalar_element = -bar_constraint_max * sig_element.numpy()
                grad_step.append(scalar_element)
            # Line search for optimizing the MLP model.
            new_parameters = LINE_SEARCH(eval_func=loss_func, para_prev=parameter_prev, grad_step=grad_step, name='Barrier loss')
            self.mlp_model.set_weights(new_parameters)

            if (np.array_equal(new_parameters, parameter_prev)):
                print("Break")
                return loss_func(new_parameters)

        return loss_func(new_parameters)
        
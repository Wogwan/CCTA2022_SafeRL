from gym.envs.registration import (
    registry,
    register,
    make,
    spec,
    load_env_plugins as _load_env_plugins,
)


register(
    id='Pygame-v0',
    entry_point='gym_game.envs:CustomEnv',
    max_episode_steps=200,
)

register(
    id="Polynomial-v0",
    entry_point="gym.envs.classic_control:PolynomialEnv",
    max_episode_steps=200,
)
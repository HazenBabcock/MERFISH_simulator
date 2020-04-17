#!/usr/bin/env python

import dotenv
import os


envPath = os.path.join(os.path.expanduser('~'), '.merlinenv')

if os.path.exists(envPath):
    dotenv.load_dotenv(envPath)

    try:
        PARAMETERS_HOME = os.path.expanduser(os.environ.get('PARAMETERS_HOME'))
        SIMULATION_PARAMETERS_HOME = os.sep.join(
                [PARAMETERS_HOME, 'simulation'])

    except TypeError:
        print('MERlin environment appears corrupt. Please run ' +
              '\'merlin --configure .\' in order to configure the environment.')
else:
    print(('Unable to find MERlin environment file at %s. Please run ' +
          '\'merlin --configure .\' in order to configure the environment.')
          % envPath)

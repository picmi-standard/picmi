import shutil
import re

# Read example script
with open('../Examples/plasma_acceleration/laser_acceleration_PICMI.py') as f:
    text = f.read()

# Replace the example code with the dummy code
new_text = re.sub('from .* import picmi', 'from dummy_code import picmi', text)

# Write the new example
with open('laser_acceleration_PICMI.py', 'w') as f:
    f.write(new_text)

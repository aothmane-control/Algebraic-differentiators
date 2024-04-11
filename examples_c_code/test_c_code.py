import os
import matplotlib.pyplot as plt

# Define the path to the files
bin_dir = ""
signal_file = os.path.join(bin_dir, "signal.txt")
result_file = os.path.join(bin_dir, "result.txt")

def read_signal_file(filename):
    """Reads the signal file and returns time and signal data."""
    time = []
    signal = []
    with open(filename, 'r') as file:
        for line in file:
            cols = line.strip().split(',')
            time.append(float(cols[0]))
            signal.append(float(cols[1]))
    return time, signal

# Read the signal files
time_signal, signal_data = read_signal_file(signal_file)
time_result, result_data = read_signal_file(result_file)

# Plot the signals
plt.figure(figsize=(10, 6))
plt.plot(time_signal, signal_data, label='Signal')
plt.plot(time_result, result_data, label='Result')
plt.xlabel('Time')
plt.ylabel('Amplitude')
plt.title('Signal vs Result')
plt.legend()
plt.grid(True)
plt.show()

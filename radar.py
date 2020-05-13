from matplotlib import pyplot as plt
import numpy as np
import math

## Radar Specifications
###########################
# Frequency of operation = 77GHz
# Max Range = 200m
# Range Resolution = 1 m
# Max Velocity = 100 m/s
###########################

# Define the target's initial position and velocity. Note : Velocity remains constant
R = 110.  # Must not exceed 200m
v = -20.  # Must be within [-70, +70] ,/s
c = 2.998e8  # Speed of light in m/s

## FMCW Waveform Generation

# Design the FMCW waveform by giving the specs of each of its parameters.
# Calculate the Bandwidth (B), Chirp Time (Tchirp) and Slope (slope) of the FMCW
# chirp using the requirements above.
range_res = 1.  # Radar range resolution
B = c / (2 * range_res)  # Chirp (sweep) bandwidth
R_max = 200.  # Radar max range
Tchirp = 5.5 * 2 * R_max / c  # Chirp (sweep) time
slope = B / Tchirp  # Chirp slope

# Operating carrier frequency of Radar
fc = 77e9  # carrier freq in Hz

# The number of chirps in one sequence. Its ideal to have 2^ value for the ease of running the FFT
# for Doppler Estimation.
Nd = 128  # No. of doppler cells OR #of sent periods % number of chirps

# The number of samples on each chirp.
Nr = 1024  # for length of time OR no. of range cells

# Timestamp for running the displacement scenario for every sample on each
# % chirp
t = np.linspace(0, Nd * Tchirp, Nr * Nd)  # total time for samples

# Creating the vectors for Tx, Rx and Mix based on the total samples input.
n_samples = len(t)
Tx = np.zeros(n_samples)  # transmitted signal
Rx = np.zeros(n_samples)  # received signal
Mix = np.zeros(n_samples)  # beat signal

# Similar vectors for range_covered and time delay.
r_t = np.zeros(n_samples)
td = np.zeros(n_samples)

## Signal generation and Moving Target simulation
# Running the radar scenario over the time.
for i in range(n_samples):
    # For each time stamp update the Range of the Target for constant velocity.

    r_t[i] = R + v * t[i]
    td[i] = 2 * r_t[i] / c

    # For each time sample we need update the transmitted and
    # received signal.

    Tx[i] = math.cos(2 * math.pi * (fc * t[i] + slope * (t[i] ** 2) / 2))
    Rx[i] = math.cos(2 * math.pi * (fc * (t[i] - td[i]) + slope * ((t[i] - td[i]) ** 2) / 2))

    # Now by mixing the Transmit and Receive generate the beat signal
    # This is done by element wise matrix multiplication of Transmit and
    # Receiver Signal

Mix = Tx * Rx

## RANGE MEASUREMENT

# reshape the vector into Nr*Nd array. Nr and Nd here would also define the size of
# Range and Doppler FFT respectively.

# Mix = np.reshape(Mix, [1, (Nr * Nd)])

# run the FFT on the beat signal along the range bins dimension (Nr) and
# normalize.

# Take the absolute value of FFT output
Mix_freq = np.abs(np.fft.fft(Mix, Nr) / Nr)

# Output of FFT is double sided signal, but we are interested in only one side of the spectrum.
# Hence we throw out half of the samples.
Mix_freq = Mix_freq[0:Nr // 2]

assert Nr % 2 == 0

# plotting the range
# figure ('Name','Range from First FFT')
# subplot(2,1,1)

# plot FFT output
plt.plot(Mix_freq)
plt.show()

## RANGE DOPPLER RESPONSE
# The 2D FFT implementation is already provided here. This will run a 2DFFT
# on the mixed signal (beat signal) output and generate a range doppler
# map.You will implement CFAR on the generated RDM

# Range Doppler Map Generation.

# The output of the 2D FFT is an image that has reponse in the range and
# doppler FFT bins. So, it is important to convert the axis from bin sizes
# to range and doppler based on their Max values.
Mix = np.reshape(Mix, [Nr, Nd], order='F')  # Note: Matlab default is to reshape by columns, Numpy default is by rows

# 2D FFT using the FFT size for both dimensions.
# sig_fft2 = np.fft.fft2(Mix,Nr,Nd);
sig_fft2 = np.fft.fft2(Mix, [Nr, Nd])

# Taking just one side of signal from Range dimension.
sig_fft2 = sig_fft2[0:Nr // 2, 0:Nd]
sig_fft2 = np.fft.fftshift(sig_fft2)
RDM = np.abs(sig_fft2)
RDM = 10 * np.log10(RDM)

# use the surf function to plot the output of 2DFFT and to show axis in both
# dimensions
doppler_axis = np.linspace(-100, 100, Nd)
range_axis = np.linspace(-200, 200, Nr // 2) * ((Nr // 2) / 400)
fig = plt.figure()
ax = fig.gca(projection='3d')

# Plot the surface.
doppler_axis, range_axis = np.meshgrid(doppler_axis, range_axis)
surf = ax.plot_surface(doppler_axis, range_axis, RDM, cmap='coolwarm',
                       linewidth=0, antialiased=True)
plt.show()

## CFAR implementation

# Slide Window through the complete Range Doppler Map

# Select the number of Training Cells in both the dimensions.

Tr = 8
Td = 4

# Select the number of Guard Cells in both dimensions around the Cell under
# test (CUT) for accurate estimation

Gr = 4
Gd = 2

# offset the threshold by SNR value in dB

offset = 5

# Create a vector to store noise_level for each iteration on training cells

n_rows = Nr // 2 - 2 * (Td + Gd)
n_cols = Nd - 2 * (Tr + Gr)
noise_level = np.zeros([n_rows, n_cols])

# design a loop such that it slides the CUT across range doppler map by
# giving margins at the edges for Training and Guard Cel ls.
# For every iteration sum the signal level within all the training
# cells. To sum convert the value from logarithmic to linear using db2pow
# function. Average the summed values for all of the training
# cells used. After averaging convert it back to logarithimic using pow2db.
# Further add the offset to it to determine the threshold. Next, compare the
# signal under CUT with this threshold. If the CUT level > threshold assign
# it a value of 1, else equate it to 0.

CFAR = np.zeros_like(RDM)
non_training_count = (2 * Gr + 1) * (2 * Gd + 1)


# Use RDM[x,y] as the matrix from the output of 2D FFT for implementing
# CFAR

# Convert decibels (dB) to power
def db2pow(xdb):
    return 10. ** (xdb / 10.)


# Convert power to decibels (dB)
def pow2db(x):
    return 10 * np.log10(x)


for i in range(n_rows):
    for j in range(n_cols):
        training_cells = db2pow(RDM[i:i + 2 * (Td + Gd) + 1, j:j + 2 * (Gr + Tr) + 1])
        training_cells[Td :-Td, Tr :-Tr] = 0
        noise_level[i, j] = pow2db(np.sum(training_cells) / (np.size(training_cells) - non_training_count))
        threshold = noise_level[i, j] + offset
        CFAR[i + Td + Gd, j + Td + Gr] = 1 if RDM[i + Td + Gd, j + Td + Gr] > threshold else 0

# The process above will generate a thresholded block, which is smaller
# than the Range Doppler Map as the CUT cannot be located at the edges of
# matrix. Hence,few cells will not be thresholded. To keep the map size same
#  set those values to 0.

# - Not applicable, as I went from edge to edge already

# display the CFAR output using the Surf function like we did for Range
# Doppler Response output.

fig = plt.figure()
ax = fig.gca(projection='3d')
# Plot the surface.
surf = ax.plot_surface(doppler_axis, range_axis, CFAR, cmap='coolwarm',
                       antialiased=True)
plt.show()

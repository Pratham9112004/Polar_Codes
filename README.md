# Polar_Codes
This project investigates Polar Codes, a groundbreaking channel coding technique that achieves Shannon’s capacity for binary-input symmetric channels. Polar Codes are at the core of modern communication systems such as 5G, making them highly relevant for research and industry applications.

📌 Key Features
Encoding/Decoding using channel polarization and transformation matrices.
BPSK Modulation and simulation over AWGN channels.
Successive Cancellation (SC) and Successive Cancellation List (SCL) decoding implemented.
LLR-based optimization to improve decoding accuracy under noisy environments.
Performance validation through BER/SNR analysis and visual plots for N = 128, 256, 512, 1024.

🔬 Topics Covered
Channel Polarization Theory & Bhattacharyya Parameter
Binary Erasure Channel (BEC)
Polar Transform and Encoding
SC and SCL Decoding Techniques
Simulation Results and Capacity Proofs
Comparison with LDPC and Turbo Codes

📊 Results
The project simulates various block lengths (N) and message lengths (K), generating:
Bit Error Rate (BER)
Block Error Rate (BLER)
Success Rates for SC and SCL Decoding
Performance plots demonstrate convergence toward Shannon’s Capacity Limit as predicted by theory.

🛠 Tools & Technologies
MATLAB for simulation and visualization
Polar code matrices constructed using Kronecker products
CRC integration in SCL decoding

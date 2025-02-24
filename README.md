# SOGRAND-C

Subject to license: "**GRAND Codebase Non-Commercial Academic Research Use License 021722.pdf**"

Non-parallelized C implementation of Soft Input Soft Output (SISO) GRAND (basic and 1-line SISO ORBGRAND) in mex format so it can be run from MATLAB.

Install: mex -O SOGRAND.c

Four MATLAB sample simulations are included:

1) sim_decoder executes OBRGAND decoding without using the SO.
2) sim_product decodes product codes with SOGRAND as the component decoder. Rows and columns are processed in MATLAB and their decoding is not parallelized,
3) sim_blkSO_acc empirically evaluates the calibration, in the Brier Score sense, of the blockwise SO.
4) sim_BLER_UER evaluates the SO control of undetected error rate.
Output is recored in results.

The following should be cited in association with results from this code.

K. R. Duffy, J. Li, and M. Medard, "Capacity-achieving guessing random additive noise decoding," IEEE Trans. Inf. Theory, vol. 65, no. 7, pp. 4023–4040, 2019.

K. R. Duffy, “Ordered reliability bits guessing random additive noise decoding," Proceedings of IEEE ICASSP, 8268–8272, 2021.

K. R. Duffy, W. An, and M. Medard, "Ordered reliability bits guessing random additive noise decoding,” IEEE Trans. Signal Process., vol. 70, pp. 4528-4542, 2022.

K. Galligan, P. Yuan, M. Médard, K. R. Duffy. "Upgrade error detection to prediction with GRAND". IEEE Globecom, 1818-1823, 2023.

P. Yuan, M. Medard, K. Galligan, K. R. Duffy. "Soft-output (SO) GRAND and Iterative Decoding to Outperform LDPC Codes". IEEE Trans. Wireless Commun., 2025.

Altough this OBRGRAND implementation is serial, the algorithm is highly parallelizable. E.g. A. Riaz, Y. Alperen, F. Ercan, W. An, J. Ngo, K. Galligan, M. Medard, K. R. Duffy, R. Yazicigil. "A sub-0.8-pJ/bit universal soft detection decoder using ORBGRAND,” IEEE J. Solid-State Circuits, 2025. Moreover, product code decoding is highly parallelizable, where all rows can be decoded in parallel and all columns can be decoded in parallel, but this is not exploited her.

For further details on GRAND, see: https://www.granddecoder.mit.edu/

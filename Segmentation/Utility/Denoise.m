function [ output_ecg ] = Denoise( input_ecg, Fs )

Fn = Fs/2;                                              % Nyquist Frequency
Wp = [1  100]/Fn;                                       % Passband (Normalised)
Ws = [0.5  110]/Fn;                                     % Stopband (Normalised)
Rp = 10;                                                % Passband Ripple (dB)
Rs = 30;                                                % Stopband Ripple (dB)
[n,Ws] = cheb2ord(Wp, Ws, Rp, Rs);                      % Chebyshev Type II Order
[b,a] = cheby2(n, Rs, Ws);                              % Transfer Function Coefficients
output_ecg = filtfilt(b, a, input_ecg);

end


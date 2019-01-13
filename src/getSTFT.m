function S = getSTFT(x,w,h,q,p)
%% Short-Time Fourier Transform
% function S = stft(x,w,h,q,p)
%
% Inputs: x - Mx1 audio
%         w - window length (default: 1024)
%         h - hop size (default: floor(w/4))
%         q - window type (default: 'hann')
%         p - padding (default: 0)
% Output: S - DxN stft matrix
%
% Note: constant signal length enforced for stft-istft iteration
%
% Johannes Traa - UIUC 2011-2013

%% check inputs
if size(x,1)  == 1;          x = x';         end
if nargin < 2 || isempty(w); w = 1024;       end
if nargin < 3 || isempty(h); h = floor(w/4); end
if nargin < 4 || isempty(q); q = 'hann';     end
if nargin < 5 || isempty(p); p = 0;          end
win = window(q,w);

%% calculate STFT
M = length(x);
N = ceil((M-w)/h)+1; % # windows
S = zeros(w+p,N);
for f=1:N-1
  S(1:w,f) = x((1:w)+(f-1)*h,1).*win;
end
x_end = x((N-1)*h+1:end,1);
S(1:w,N) = [x_end; zeros(w-length(x_end),1)];
S = fft(S,w+p,1);
S = S(1:(w+p)/2+1,:);
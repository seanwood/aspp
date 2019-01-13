function x = getISTFT(S,w,h,q)
%% Inverse Short-Time Fourier Transform
% function x = istft(S,w,h,q)
% 
% Inputs: S - DxN complex-valued STFT matrix
%         w - window size (default: 1024)
%         h - hop size (default: floor(w/4))
%         q - window type (default: 'hann')
% Output: x - Lx1 audio
%
% See Griffin and Lim (1984)
%
% Johannes Traa - UIUC 2011-2013

N = size(S,2); % # windows

%% check inputs
if size(S,1) ~= w; S = [S; conj(S(end-1:-1:2,:))]; end
if nargin < 2 || isempty(w); w = 1024;       end
if nargin < 3 || isempty(h); h = floor(w/4); end
if nargin < 4 || isempty(q); q = 'hann';     end
win = window(q,w);

%% invert and overlap-add frames
L = (N-1)*h + w; % rec length
x = zeros(L,1); % rec signal
s = zeros(L,1); % sum of squared windows
y = ifft(S,[],1);
for f=1:N
  g = (f-1)*h;
  x((1:w)+g) = x((1:w)+g) + real(y(:,f)).*win;
  s((1:w)+g) = s((1:w)+g) + win.^2;
end
x(h:L-h+1) = x(h:L-h+1)./s(h:L-h+1);
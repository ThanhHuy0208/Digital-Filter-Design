                    % Đề tài 3: Thiết kế bộ lọc số
                    % GVHD    : PGS.TS Hà Hoàng Kha

% Công thức, hàm và các thuật toán được tham khảo từ: 
% [1] S. J. Orfanidis, Introduction to Signal Processing, Prentice –Hall Publisher 2010.
% [2] V. K. Ingle, J. Proakis, Digital Signal Processing Using Matlab, Cengage Learning, 3 Edt, 2011.
                    fprintf('Thiết kế bộ lọc số \n');
                    boloc = input('1: FIR, 2: IIR \n');
                    switch boloc 
    
                        case 1 % FIR

% Thiết kế bộ lọc FIR bằng phương pháp cửa sổ 
fprintf('Thiết kế bộ lọc FIR bằng phương pháp cửa sổ\n');
fprintf('Chọn bộ lọc : \n');
type = input('1: Thông thấp, 2: Thông cao, 3: Thông dải, 4: Chắn dải \n');
fprintf('Nhập các thông số cho bộ lọc \n');
Fs = input('Chọn tần số lấy mẫu (Hz): ');    % tần số lấy mẫu = Fs
delta = input('Độ gợn sóng: delta = ');      % max ripple 
% Window method -> ripple stopband = ripple passband
A = 20 * log10(delta);  % Chuyển độ gợn sóng (suy hao) sang thang dB 
fprintf('Độ gợn sóng (dB) A = %f\n', A);
%-------------------------------------------------
fprintf('Hãy nhập các tần số chuẩn hóa (x pi rad/sample) \n');
if (type == 1) || (type == 2)
wp = pi * input('Nhập tần số cắt dải thông: wp = ');
ws = pi * input('Nhập tần số cắt dải chắn : ws = ');
wc = (ws + wp)/2; 
D  = abs(ws - wp); % bề rộng dải chuyển tiếp 
elseif (type == 3) || (type == 4)
wp1 = pi * input('Nhập tần số cắt dưới dải thông: wp1 = '); 
wp2 = pi * input('Nhập tần số cắt trên dải thông: wp2 = ');
ws1 = pi * input('Nhập tần số cắt dưới dải chắn : ws1 = ');
ws2 = pi * input('Nhập tần số cắt trên dải chắn : ws2 = ');
wc1 = (ws1 + wp1)/2; % tần số cắt dưới
wc2 = (ws2 + wp2)/2; % tần số cắt trên
D   = min(abs(ws1-wp1), abs(ws2-wp2)); %  bề rộng dải chuyển tiếp
end
fprintf('Bề rộng dải chuyển tiếp D = %f\n',D);
%--------------------------------------------
% Dựa vào độ gợn sóng, chọn cửa sổ tương ứng 
% L là chiều dài bộ lọc, L-1 là bậc của bộ lọc 
% Do đối xứng chẵn nên L luôn phải là số lẻ
% Hàm window đặc trưng cho từng cửa sổ 
% Chọn hàm cửa sổ phù hợp 
if (A >= -21)
    fprintf('Chọn cửa sổ Rectangular\n'); 
    L = ceil(1.8*pi/D) + 1;
    L = L + (rem(L,2) == 0);
    window = boxcar(L);
    
elseif (A < -21)&&(A >= -44)
    fprintf('Chọn cửa sổ Hann\n'); 
    L = ceil(6.2*pi/D) + 1;
    L = L + (rem(L,2) == 0);
    window = hanning(L);
    
elseif (A < -44)&&(A >= -53)
    fprintf('Chọn cửa sổ Hamming\n');
    L = ceil(6.6*pi/D) + 1;
    L = L + (rem(L,2) == 0);
    window = hamming(L);
    
elseif (A < -53)&&( A >= -74)
    fprintf('Chọn cửa sổ Blackman\n');
    L = ceil(11*pi/D) + 1;
    L = L + (rem(L,2) == 0);
    window = blackman(L);
end 
% L = L + (rem(L,2) == 0)    % rem(L,2) là phép chia lấy số dư 
% L lẻ thì   L <-- L
% L chẵn thì L <-- L + 1  
fprintf('Chiều dài của bộ lọc: L = %d\n',L);

%----------Tạo hàm xung dirac----------------  
    n = 0 : 1 : (L-1); 
    m = n - (L-1)/2 + eps ; 
    dirac = (m == eps) ; 
%--------------------------------------------
% Chọn hàm đáp ứng xung lí tưởng dựa theo kiểu bộ lọc
% hàm ideal_lp được để dưới cùng của bài code
switch type 
    case 1 % thông thấp 
        hd = ideal_lp(wc, L); 
    case 2 % thông cao
        hd = dirac - ideal_lp(wc, L);
    case 3 % thông dải 
        hd = ideal_lp(wc2, L) - ideal_lp(wc1, L);
    case 4 % chắn dải 
        hd = dirac - (ideal_lp(wc2, L) - ideal_lp(wc1, L));
end
window = window' ; % đảo ma trận cột thành hàng 
h = hd .* window ; % nhân trong miền thời gian 
%-----Bây giờ h mới là đáp ứng xung thực tế của bộ lọc---------
% Vẽ các đồ thị 
% Đáp ứng xung lí tưởng 
subplot(1,1,1);
subplot(2,2,1); stem(n,hd); title('Đáp ứng xung lí tưởng'); 
axis([0 L-1 -0.2 0.5]); xlabel('n'); ylabel('hd(n)');
% Hàm cửa sổ đã chọn 
subplot(2,2,2); stem(n, window); title('Hàm cửa sổ đã dùng'); 
axis([0 L-1 0 1.1]); xlabel ('n'); ylabel ('w(n)');
% Đáp ứng xung thực tế
subplot(2,2,3); stem(n, h); title('Đáp ứng xung thực tế'); 
axis([0 L-1 -0.2 0.5]); xlabel('n'); ylabel('h(n)');
% Đáp ứng biên độ (dB)
[H,w] = freqz(h, 1);
mag = abs(H);                         % chuyển sang độ lớn biên độ 
Hdb = 20*log10( (mag+eps)/max(mag) );   % chuyển sang dB (+eps tránh log0)
subplot(2,2,4); plot(w/pi, Hdb); title('Đáp ứng biên độ ở thang dB'); grid on;
axis([0 1 -100 10]); 
xlabel('Tần số chuẩn hóa (\times\pi rad/sample)'); ylabel('Biên độ (dB)');
%Đáp ứng biên độ (Linear Scale) 
figure; 
plot(w/pi, mag); title('Đáp ứng biên độ ở thang tuyến tính (Linear Scale)'); grid on;
axis([0 1 -0.3 1.3]); 
xlabel('Tần số chuẩn hóa (\times\pi rad/sample)'); ylabel('Biên độ ');

                        case 2 % IIR

% Thiết kế bộ lọc số IIR bằng phương pháp bilinear transformation.
fprintf('Thiết kế bộ lọc IIR bằng phương pháp bilinear transformation \n');
fprintf('Chọn bộ lọc : \n');
type = input('1: Thông thấp, 2: Thông cao, 3: Thông dải, 4: Chắn dải \n');
% Chọn loại bộ lọc
tp = input('1: Butterworth, 2: Chebyshev 1 , 3: Chebyshev 2, 4: Elliptic \n');
fprintf('Nhập các thông số cho bộ lọc \n');
Fs = input('Chọn tần số lấy mẫu (Hz): '); 
Rp = input('Độ gợn sóng dải thông (dB): Rp = '); % ripple passband
    epsilon_p = sqrt((10^(Rp/10))-1) ; 
    fprintf('epsilon_p = %f\n', epsilon_p);
As = input('Độ suy hao dải chắn   (dB): As = '); % attenuation stopband
    epsilon_s = sqrt((10^(As/10))-1) ;
    fprintf('epsilon_s = %f\n', epsilon_s);
%-------------------------------------------------
% Nhập các tần số miền số theo yêu cầu đề bài, 
% từ đó chuyển các thông số qua miền tương tự để tính toán
fprintf('Hãy nhập các tần số chuẩn hóa (x pi rad/sample) \n');
if (type == 1) || (type == 2)
wp = pi * input('Nhập tần số cắt dải thông: wp = ');
ws = pi * input('Nhập tần số cắt dải chắn : ws = ');
    if   (type == 1)
        % Chuyển qua miền tương tự
        omega_p = tan(0.5*wp) ; 
        omega_s = tan(0.5*ws) ;
    else (type == 2)
        % Chuyển qua miền tương tự
        omega_p = (1 / tan(0.5*wp)) ; 
        omega_s = (1 / tan(0.5*ws)) ;
    end 
elseif (type == 3) || (type == 4)
wp1 = pi * input('Nhập tần số cắt dưới dải thông: wp1 =  '); 
wp2 = pi * input('Nhập tần số cắt trên dải thông: wp2 =  ');
ws1 = pi * input('Nhập tần số cắt dưới dải chắn : ws1 =  ');
ws2 = pi * input('Nhập tần số cắt trên dải chắn : ws2 =  ');
wp=[wp1 wp2];
ws=[ws1 ws2];
c = sin(wp1 + wp2) / (sin(wp1) + sin(wp2)) ;   
    if (type == 4)
         % Chuyển qua miền tương tự
        omega_s1 = sin(ws1) / (cos(ws1) - c) ;
        omega_s2 = sin(ws2) / (cos(ws2) - c) ;
        omega_p  = abs( sin(wp2) / (cos(wp2) - c) ) ; 
        
    else (type == 3)
         % Chuyển qua miền tương tự
        omega_s1 = (c - cos(ws1)) / sin(ws1) ;
        omega_s2 = (c - cos(ws2)) / sin(ws2) ;
        omega_p  = abs( (c - cos(wp2)) / sin(wp2) ) ; 
    end 
omega_s  = min(abs(omega_s1) , abs(omega_s2));
end
%--------------------------------------------------------------------------
% Bậc của bộ lọc 
    N = ceil( (log(epsilon_s/epsilon_p)) / (log(omega_s/omega_p)) );
    % Tần số cắt -3dB trong miền tương tự
    omega_c = omega_p / (epsilon_p ^ (1/N) ) ;        
    fprintf('Sau tính toán, bậc của bộ lọc là %d\n', N); 
%--------------------------------------------------------------------------
% Tìm tần số cắt -3dB của bộ lọc trong miền số 
if (type == 1)  % Bộ lọc thông thấp 
    % Tần số cắt -3dB trong miền số (Hz)
    fc = (Fs/pi)* atan(omega_c) ;
    % Tần số cắt -3dB trong miền số (radian/sample)
    wn = (2*pi*(fc/Fs))/pi;
end
if (type == 2)  % Bộ lọc thông cao 
    % Tần số cắt -3dB trong miền số (Hz)    
    fc = (Fs/pi)* atan(1/omega_c) ; 
    % Tần số cắt -3dB trong miền số (radian/sample)
    wn = (2*pi*(fc/Fs))/pi;
end
if (type == 3) % Bộ lọc thông dải 
    w0a= 2*(atan( (sqrt(omega_c^2 + 1 - c^2) - omega_c) / (1 + c)));
    w0b= 2*(atan( (sqrt(omega_c^2 + 1 - c^2) + omega_c) / (1 + c)));
    % Tần số cắt thấp trong miền số (Hz)
    fca=(Fs*w0a)/(2*pi) ; 
    % Tần số cắt cao trong miền số (Hz)
    fcb=(Fs*w0b)/(2*pi) ; 
    % Tần số cắt thấp trong miền số (radian/sample)
    wna = w0a/pi;
    % Tần số cắt cao trong miền số (radian/sample)
    wnb = w0b/pi;
    wn = [wna wnb];
end
if (type == 4) % Bộ lọc chắn dải
    w0a = 2*(atan((sqrt( 1+ (omega_c^2)*(1- c^2))-1) /(omega_c* (1 + c))));
    w0b = 2*(atan((sqrt( 1+ (omega_c^2)*(1- c^2))+1) /(omega_c* (1 + c))));
    % Tần số cắt thấp trong miền số (Hz)
    fca=(Fs*w0a)/(2*pi) ; 
    % Tần số cắt cao trong miền số (Hz)
    fcb=(Fs*w0b)/(2*pi) ; 
    % Tần số cắt thấp trong miền số (radian/sample)
    wna = w0a/pi;
    % Tần số cắt cao trong miền số (radian/sample)
    wnb = w0b/pi;
    wn = [wna wnb];
end
%----------------------------------------------------------------------
% Thiết kế bộ lọc khi đã có bậc bộ lọc N và tần số cắt wn (rad/sample)
if(tp==1)   % Bộ lọc Butterworth
switch(type)
    case 1
         [z,p,k]= butter(N,wn,"low");
    case 2
         [z,p,k]= butter(N,wn,"high");
    case 3
         [z,p,k]= butter(N,wn,"bandpass");
    case 4 
         [z,p,k]= butter(N,wn,"stop");
end
end
if(tp==2)    % Bộ lọc Chebyshev loại I
switch(type)
    case 1
         [z,p,k]= cheby1(N,Rp,wn,"low");
    case 2
         [z,p,k]= cheby1(N,Rp,wn,"high");
    case 3
         [z,p,k]= cheby1(N,Rp,wn,"bandpass");
    case 4 
         [z,p,k]= cheby1(N,Rp,wn,"stop");
end
end
if(tp==3)    % Bộ lọc Chebyshev loại II
switch(type)
    case 1
         [z,p,k]= cheby2(N,As,wn,"low");
    case 2
         [z,p,k]= cheby2(N,As,wn,"high");
    case 3
         [z,p,k]= cheby2(N,As,wn,"bandpass");
    case 4 
         [z,p,k]= cheby2(N,As,wn,"stop");
end
end
if(tp==4)    % Bộ lọc Elliptic 
switch(type)
    case 1
        [z,p,k]= ellip(N,Rp,As,wn,"low");
    case 2
        [z,p,k]= ellip(N,Rp,As,wn,"high");
    case 3
        [z,p,k]= ellip(N,Rp,As,wn,"bandpass");
    case 4 
        [z,p,k]= ellip(N,Rp,As,wn,"stop");
end
end
%--------------------------------------------------------------------------
% Vẽ đáp ứng biên độ và pha ở thang dB
    sos = zp2sos(z,p,k);
    freqz(sos);
% Vẽ đáp ứng biên độ ở thang tuyến tính
    [H,w] = freqz(sos);
    mag = abs(H);
%Đáp ứng biên độ (Linear Scale) 
figure; 
plot(w/pi, mag); title('Đáp ứng biên độ ở thang tuyến tính (Linear Scale)'); grid on;
axis([0 1 -0.3 1.3]); 
xlabel('Tần số chuẩn hóa (\times\pi rad/sample)'); ylabel('Biên độ ');

                           end 

%----------------------------------------------------------------------------
%-----------------Khảo sát các bộ lọc bằng cách đưa các tín hiệu ngõ vào --------
%----------------------------------------------------------------------------
% Chuyển sang tần số f(Hz) để khảo sát
% In ra màn hình các tần số cắt f(Hz) của bộ lọc vừa thiết kế
if (type == 1) || (type == 2)
        fp = (wp * Fs) / (2*pi);
        fs = (ws * Fs) / (2*pi);
        fprintf('Tần số cắt dải chắn của bộ lọc:  fs(Hz) = %f\n',fs);
        fprintf('Tần số cắt dải thông của bộ lọc: fp(Hz) = %f\n',fp);
elseif (type == 3) || (type == 4)
        fp1 = (wp1 * Fs) / (2*pi);
        fs1 = (ws1 * Fs) / (2*pi);
        fp2 = (wp2 * Fs) / (2*pi);
        fs2 = (ws2 * Fs) / (2*pi);
    fprintf('Tần số cắt dưới dải chắn của bộ lọc : fs1(Hz) = %f\n',fs1);
    fprintf('Tần số cắt dưới dải thông của bộ lọc: fp1(Hz) = %f\n',fp1);
    fprintf('Tần số cắt trên dải chắn của bộ lọc : fs2(Hz) = %f\n',fs2);
    fprintf('Tần số cắt trên dải thông của bộ lọc: fp2(Hz) = %f\n',fp2);
end
% Nhập số lượng tín hiệu đầu vào cần lọc.
signals = input('Nhập số lượng tín hiệu đầu vào cần lọc : ');

% Tạo mảng chứa các tần số và biên độ của các tín hiệu sắp nhập.
fin  = zeros(1, signals);   % chứa các tần số f(Hz) đi vào bộ lọc
ain  = zeros(1, signals);   % chứa các biên độ đi vào bộ lọc
win  = zeros(1, signals);   % Tần số f(Hz) chuyển thành w (rad/sample)
aout = zeros(1, signals);   % chứa các biên độ sau khi lọc

% Chọn các thông số cho từng tín hiệu.
for i = 1:signals
    fprintf('Tín hiệu thứ %d:\n', i);
    fin(i)  = input('  Nhập tần số f(Hz) = ');
    ain(i)  = input('  Nhập biên độ A = ');
    win(i)  = (2*fin(i))/Fs; % chuyển tần số f(Hz) vừa nhập thành w để xử lí
    aout(i) = ain(i)*mag(round(win(i)*512)); % Hàm freqz chia thành 512 khoảng
    % Lấy biên độ vừa nhập nhân với biên độ của đáp ứng xung ngay tại tần
    % số mong muốn --> ta được biên độ của tín hiệu sau khi lọc
end
% Vẽ phổ dương của các tín hiệu trước khi lọc
figure; 
stem(fin, ain, 'filled');
title('Phổ tín hiệu trước khi lọc');
axis([min(fin)-1000 max(fin)+1000 -1 max(ain)+1]);
xlabel('Tần số (Hz)'); ylabel('Biên độ');
grid on;
% Vẽ phổ dương của các tín hiệu sau khi lọc
figure; 
stem(fin, aout, 'filled'); 
title('Phổ tín hiệu sau khi lọc');
axis([min(fin)-1000 max(fin)+1000 -1 max(ain)+1]);
xlabel('Tần số (Hz)'); ylabel('Biên độ');
grid on; 

% Hàm đáp ứng xung dưới đây được dùng trong bài FIR
%-----Hàm đáp ứng xung lí tưởng  
function hd = ideal_lp(wc, L) % hd(n) cho LPF
    alpha = (L-1)/2 ;           % độ dời 
    n = 0 : 1 : (L-1);
    % Để tránh n = alpha, cộng thêm 1 số rất nhỏ (eps) 
    m = n - alpha + eps ;       % tránh chia cho 0
    hd = sin(wc*m)./ (pi*m);    % đáp ứng xung lí tưởng của LPF
  % hd = (m == 0)*(wc/pi) + (m != 0)*(sin(wc*m)./(pi*m)) ;
end



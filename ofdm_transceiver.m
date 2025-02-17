clc; clear; 
%% PARAMETERS
ch_bandwidth = 10e6;
n_fft = 1024; 
no_occupied_sc = 600;
sc_spacing = 15e3;
no_rb = 50;
code_rate = 1/3;
mod_order = 64; % 64 QAM
bits_per_symbol = log2(mod_order); % 6
first_cp_length = 80;
cp_length = 72;
pilot_sc_spacing = 20;
no_pilots_sc = no_occupied_sc/pilot_sc_spacing;

SNR_dB = 0:5:30;
BER = zeros(size(SNR_dB));
tx_codeword = [];

img = imread('woman_darkhair.tif');  % Load the image
[rows, columns, numberOfColorChannels] = size(img);
if numberOfColorChannels > 1
    img_gray = rgb2gray(img);  % Convert to grayscale
else
    img_gray = img;
end
img_bits = de2bi(img_gray(:), 8, 'left-msb'); % Convert pixel values to 8-bit binary
img_bits = img_bits(:); 
data_size = length(img_bits); 
figure
tiledlayout(2,4);

nexttile
imshow(img_gray)
title("Original image")

%% TRANSMITTER
% Perform Turbo Encoding
n = 1;
K = 6144;
intrlvrIndices = helperLTEIntrlvrIndices(K);
turboEnc = comm.TurboEncoder('TrellisStructure', poly2trellis(4, [13 15], 13),'InterleaverIndices', intrlvrIndices);

while (n+6143) < data_size
    data = img_bits(n:n+6143);
    data_coded = turboEnc(data);
    tx_codeword = [tx_codeword; data_coded];
    n = n + 6144;
end

data = img_bits(n:data_size);
intrlvrIndices = helperLTEIntrlvrIndices(length(data));
turboEnc = comm.TurboEncoder('TrellisStructure', poly2trellis(4, [13 15], 13),'InterleaverIndices', intrlvrIndices);
data_coded = turboEnc(data);
tx_codeword = [tx_codeword; data_coded];

% Symbol Mapping and Modulation
tx_symbols = lteSymbolModulate(tx_codeword,'64QAM');

% Pilot insertion
no_data = length(tx_symbols);
no_pilots = floor(no_data / (pilot_sc_spacing - 1)); 
num_total_symbols = no_data + no_pilots;
tx_total_bits = zeros(num_total_symbols, 1); 
pilot_positions = 20:20:num_total_symbols;  
pilot_positions(pilot_positions > num_total_symbols) = [];  
pilot = 1 + 1i;
tx_total_bits(pilot_positions) = pilot;
data_positions = setdiff(1:num_total_symbols, pilot_positions);
tx_total_bits(data_positions) = tx_symbols;

% TX Resource grid
% padding 
if rem(length(tx_total_bits), no_occupied_sc) ~= 0
    padding_length = no_occupied_sc - rem(length(tx_total_bits), no_occupied_sc);
    padding = zeros(padding_length, 1) + 1i * 0; 
    tx_total_bits = [tx_total_bits; padding];
end

num_ofdm_symbols = ceil(length(tx_total_bits)/no_occupied_sc);
symbol_with_dc = zeros(no_occupied_sc + 1, 1);
tx_resource_grid = zeros(n_fft,num_ofdm_symbols);
column_index = 1;
for i = 1:600:length(tx_total_bits)
    symbol_bits = tx_total_bits(i:i + no_occupied_sc - 1);
    dc_position = 301;
    symbol_with_dc = [symbol_bits(1:300); 0; symbol_bits(301:end)];
    tx_resource_grid(214:214 + length(symbol_with_dc) - 1, column_index) = symbol_with_dc;
    column_index = column_index + 1;
end

% IFFT
ofdm_symbols_ifft = ifft(tx_resource_grid,n_fft);
ofdm_symbols = ofdm_symbols_ifft(:);


% CP addition
tx_ofdm_with_cp = zeros(((n_fft + first_cp_length) + (num_ofdm_symbols-1) * (n_fft + cp_length)),1); % 2013360 x 1
k = 1;
for i = 1:num_ofdm_symbols
    symbol = ofdm_symbols((i-1)*n_fft + 1:i*n_fft);
    if i == 1
        cp = symbol(end-first_cp_length+1:end);
        symbol_with_cp = [cp; symbol];
        tx_ofdm_with_cp(k:k + n_fft + first_cp_length - 1) = symbol_with_cp;
        k = k + n_fft + first_cp_length;
    else
        cp = symbol(end-cp_length+1:end);
        symbol_with_cp = [cp; symbol];
        tx_ofdm_with_cp(k:k + n_fft + cp_length - 1) = symbol_with_cp;
        k = k + n_fft + cp_length;
    end
end

%% CHANNEL AND RECEIVER
for snr_idx = 1: length(SNR_dB)
    snr = SNR_dB(snr_idx);

    % Channel
    signalPower = mean(abs(tx_ofdm_with_cp).^2);
    rx_ofdm_signal = awgn(tx_ofdm_with_cp,snr, 'measured');

    % Cp removal and FFT
    rx_resource_grid = zeros(n_fft, num_ofdm_symbols);
    for i = 1:num_ofdm_symbols
        if i == 1     
            symbol = rx_ofdm_signal(first_cp_length+1:first_cp_length + 1024); 
        else
            symbol = rx_ofdm_signal((i-1) * (n_fft + cp_length)+ first_cp_length + 1:(i-1) * (n_fft + cp_length)+ first_cp_length + 1024);
        end
        rx_resource_grid(:, i) = fft(symbol, n_fft);
    end

    % Channel estimation using LS
    rx_symbols = [];
    for i = 1:num_ofdm_symbols
        if i == num_ofdm_symbols
            symbol = rx_resource_grid(214:214 + no_occupied_sc - length(padding),i);
        else
            symbol = rx_resource_grid(214:214 + no_occupied_sc,i);
        end
        symbol(dc_position) = []; % removing the dc component
        active_subcarriers = (1:length(symbol))';
        pilot_positions = 20:pilot_sc_spacing:length(symbol);
        received_pilots = symbol(pilot_positions);
        H_est_pilots = received_pilots ./ pilot;
        H_est = interp1(pilot_positions, H_est_pilots, active_subcarriers, 'linear', 'extrap');
        rx_data_symbols = symbol ./ H_est;
        % removing the pilots
        data_positions = setdiff(1:length(symbol), pilot_positions);
        rx_symbols_cur = rx_data_symbols(data_positions);
        rx_symbols = [rx_symbols;rx_symbols_cur];
    end
    
    % Demodulation
    rx_codeword = lteSymbolDemodulate(rx_symbols,'64QAM','Soft');

    % Turbo Decoding
    rx_bits = [];
    K = 6144;
    block_size = 3*K + 12;
    intrlvrIndices = helperLTEIntrlvrIndices(K);
    turboDec  = comm.TurboDecoder('TrellisStructure', poly2trellis(4, [13 15], 13),'InterleaverIndices', intrlvrIndices);
    n = 1;
    while (n + block_size - 1) < length(rx_codeword)
        data = rx_codeword(n:n + block_size - 1);
        data_decoded = turboDec(data);
        rx_bits = [rx_bits;data_decoded];
        n = n + block_size;
    end
    
    K = 2048;
    data = rx_codeword(n:length(rx_codeword));
    intrlvrIndices = helperLTEIntrlvrIndices(K);
    turboDec  = comm.TurboDecoder('TrellisStructure', poly2trellis(4, [13 15], 13),'InterleaverIndices', intrlvrIndices);
    data_decoded = turboDec(data);
    rx_bits = [rx_bits;data_decoded]; % 2097152 x 1
    
    % BER calculation
    errors = sum(rx_bits ~= img_bits);
    BER(snr_idx) = errors / length(img_bits);
    fprintf("SNR = %d dB, BER = %e\n", snr, BER(snr_idx));
    
    % Image reconstruction
    rx_bits_matrix = reshape(rx_bits, [], 8); 
    rx_pixel_values = bi2de(rx_bits_matrix, 'left-msb');
    rx_image = reshape(rx_pixel_values, [512, 512]);
    
    nexttile
    imshow(uint8(rx_image))
    title(sprintf('Recived image at SNR = %d dB', SNR_dB(snr_idx)));
end

figure;
semilogy(SNR_dB, BER, 'bo-', 'LineWidth', 1, 'MarkerSize', 4);
grid on;
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
title('BER vs SNR for OFDM with Turbo Coding');
legend('64-QAM, Code Rate 1/3');

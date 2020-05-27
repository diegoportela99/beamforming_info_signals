% MISO Beamforming simulation

%initial parameters
symbol_count = (10)^(6);
modulation_order = 2;
modulation_points = (2)^(modulation_order);
signal_to_noise = (-25):(1):(35);
symbol_ETO = (signal_to_noise + 3)*(modulation_order-1); %Energy-to-Noise
linear_SNR = 10.^(signal_to_noise / 10);
estimation_error = zeros(1, length(signal_to_noise)); %vector set to zero
intended_receiver_error = estimation_error;
nonintended_receiver_error = estimation_error;
timer_start = tic;

for loop = (1): (1):(length(signal_to_noise))
    % develop randomised channel coefficents
    tx_to_intended_channel = repelem(reshape((randn(1,symbol_count) + ...
        randn(1,symbol_count)*1j) / sqrt(2),2,symbol_count/2),1,2);

    tx_to_nonintended_channel = repelem(reshape((randn(1,symbol_count) ...
        + randn(1,symbol_count)*1j) / sqrt(2),2,symbol_count/2),1,2);
    
    % nonintended transmit QPSK Grey-coded modulation
    x = round(rand(1,symbol_count)) + round(rand(1,symbol_count)) * 2;
    b = reshape(dec2bin(x).',1,2 * symbol_count);
    
    %init store transmit symbol vector
    store = zeros(1,symbol_count); %vector with symbol size set to zero
    for size_u = (1 : symbol_count)
        switch(x(size_u))
            case 0
                store(size_u) = -1;
            case 1
                store(size_u) = -1j;
            case 2
                store(size_u) = 1j;
             ...
            otherwise
                store(size_u) = 1;
        end
    end
 
    %creating noise in channels
    tx_to_intended_noise = 10^(-symbol_ETO(loop)/20) * ...
    (randn(1,symbol_count) + randn(1,symbol_count)*1j)/sqrt(2);

    tx_to_nonintended_noise = 10^(-symbol_ETO(loop)/20) * ...
    (randn(1,symbol_count) + randn(1,symbol_count)*1j)/sqrt(2);

    %beamforming
    tx_to_intended_beam = tx_to_intended_channel.*exp(-1j * ...
        angle(tx_to_intended_channel));
    
    %transmission
    store = repelem(store, 2, 1) / sqrt(2); %set tx symbols
    
    intended_reciever = sum(tx_to_intended_beam.*store, 1) + ...
        tx_to_intended_noise;
    nonintended_reciever = sum(tx_to_nonintended_channel.*store, 1) + ...
        tx_to_nonintended_noise;
    
    %Equalization for intended channel - beamforming
    intended_s_estimate = intended_reciever./sum(tx_to_intended_beam,1);
    
    %symb detect - intended with beamforming to intended
    angle_intended = angle(intended_s_estimate) * 180/pi;
    [int_x_est, int_b_est] = symb_detect(angle_intended);
    
    %nonintended with beamforming to intended
    nonintended_s_estimate = nonintended_reciever./ ...
        sum(tx_to_nonintended_channel,1);
    
    %symb detect - nonintended with Beamforming to intended
    angle_nonintended = 180/pi * angle(nonintended_s_estimate);
    [non_x_est, non_b_est] = symb_detect(angle_nonintended);
    
    % Count the estimation errors
    intended_receiver_error(loop) = size(find(b - int_b_est),2);
    nonintended_receiver_error(loop) = size(find(b - non_b_est),2);
    
    % Display elapsed time
    display(toc(timer_start));
end

% Simulation results
intended_BER = intended_receiver_error/(2 * symbol_count);
nonintended_BER = nonintended_receiver_error/(2 * symbol_count);

% Plot the results
close all
figure
hold on
semilogy(signal_to_noise,intended_BER,'.-c','LineWidth',1); % MISO (intended - BF)
semilogy(signal_to_noise,nonintended_BER,'-ro','LineWidth',1); % MISO (non intended - BF)
axis([-25 35 10^-5 1])
grid on
title('QPSK Bit Error Rate (BER) - Beamforming');
legend('intended (nTx=2, nRx=1, Tx BF)','non intended (nTx=2, nRx=1, Tx BF)');
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');


%detection of symbol - function
function [x_estimate, b_estimate] = symb_detect(angle)
    symbol_count = (10)^(6);
    x_estimate = zeros(1,symbol_count);
    for size_i = (1 : symbol_count)
         if -45 <= angle(size_i) && angle(size_i) < 45
             x_estimate(size_i) = 3;
             
         elseif 45 <= angle(size_i) && angle(size_i) < 135
             x_estimate(size_i) = 2;
             
         elseif 135 <= angle(size_i) || -135 >= angle(size_i)
             x_estimate(size_i) = 0;
             
         else
             x_estimate(size_i) = 1;
         end
    end
    b_estimate = reshape(dec2bin(x_estimate).',1,2 ...
        * symbol_count);
end


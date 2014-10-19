for idx = 1:3
    
    % Generate strings for file and image names:
    file_name = sprintf('inflammation-%02d.csv', idx);
    img_name = sprintf ('patient_data-%02d.png', idx);

    patient_data = csvread(file_name);
    ave_inflammation = mean(patient_data, 1);

    figure()

    subplot(2, 2, 1);
    plot(ave_inflammation);
    ylabel('average')

    subplot(2, 2, 2);
    plot(max(patient_data, [], 1));
    ylabel('max')

    subplot(2, 2, 3);
    plot(min(patient_data, [], 1));
    ylabel('min')

    print('-dpng', img_name);
    close();

end

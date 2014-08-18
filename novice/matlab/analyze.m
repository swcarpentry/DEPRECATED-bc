for i = 1:3
    
    % Generate strings for file and image names:
    file_name = sprintf('inflammation-%02d.csv', i);
    img_name = sprintf ('patient_data-%02d.svg', i);

    patient_data = csvread(file_name);
    ave_inflammation = mean(patient_data, 1);

    figure()

    subplot(1, 3, 1);
    plot(ave_inflammation);
    ylabel('average')

    subplot(1, 3, 2);
    plot(max(patient_data, [], 1));
    ylabel('max')

    subplot(1, 3, 3);
    plot(min(patient_data, [], 1));
    ylabel('min')

    print(img_name);
    close();

end

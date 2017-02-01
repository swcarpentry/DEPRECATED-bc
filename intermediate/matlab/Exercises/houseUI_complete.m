function houseUI_complete()
    
    % Load data.
    handles = load('Ex04_House.mat');
    
    % Create background figure.
    screenPos = get(0, 'ScreenSize');
    dxScr = screenPos(3); dyScr = screenPos(4);
    handles.fh = figure('Name', 'House Price to Median Earnings Ratio UI', ...
                        'Position', [dxScr/4, dyScr/4, dxScr/2, dyScr/2], ...
                        'NumberTitle', 'off', ...
                        'Color', 0.9*ones(1, 3));
    
    % Create axes.
    handles.ah = axes('Parent', handles.fh, ...
                      'Position', [0.5, 0.3, 0.4, 0.6], ...
                      'XLim', [1997, 2012]);                      
    
    % Create pop-up menu.
    handles.popuph = uicontrol(handles.fh, ...
        'Style', 'popupmenu', ...
        'Units', 'normalized', ...
        'Position', [0.05, 0.1, 0.3, 0.6], ...
        'String', handles.pricesToWages.Properties.RowNames, ...
        'BackgroundColor', 0.85*ones(1, 3), ...
        'Callback', @plotPrices, ...
        'TooltipString', 'Choose region of interest from this list');
        
    function plotPrices(varargin)
        % TODO: Complete the code here so that it updates the plot when a 
        % selection is made from the pop-up menu.
        
        % Required row number.
        k = get(handles.popuph, 'Value');
        
        % Update the plot with the correct data.
        plot(handles.ah, 1997:2012, handles.pricesToWages{k, :}, ...
            'LineWidth', 2, 'Marker', '*')
        xlabel(handles.ah, 'Year')
        ylabel(handles.ah, 'House Price to Median Earnings Ratio')
        title(handles.ah, ...
            sprintf('Data for %s', handles.pricesToWages.Properties.RowNames{k}), ...
            'FontWeight', 'Bold')
        grid(handles.ah)

        
    end % plotPrices
    
    
end % houseUI_complete
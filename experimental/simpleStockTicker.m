classdef simpleStockTicker < handle
% simpleStockTicker   A simple stock ticker app written with a MATLAB class
%
% This app provides an example of how to use a MATLAB class to write a
% simple stock ticker which updates a plot of stock prices over time for a
% given ticker symbol.
%
% David Garrison @ The MathWorks


    properties
        Figure                  % Graphics handles
        Axis
        Line
        TickerText
        TickerEdit

        Timer                   % Timer object to get updated prices
        TimerUpdateRate = 1     % In seconds
        NumValues = 30          % Number of values shown in the plot
        TickerSymbol = 'GOOG'   % Current ticker symbol (initial value = 'GOOG')
    end

    methods

        function app = simpleStockTicker
        % This is the "constructor" for the class
        % It runs when an object of this class is created
            app.Figure = figure('MenuBar','none',...           % Main figure
                'NumberTitle','off','Name','Simple Stock Ticker',...
                'CloseRequestFcn',@app.closeApp) ;
            app.Axis = axes('Parent',app.Figure,...            % Axis for prices
                'Position',[.13 .15 .78 .75]);
            app.TickerText = uicontrol(app.Figure,...          % 'Symbol' label
                'Style','text','Position',[20 20 50 20],...
                'String','Symbol:');
            app.TickerEdit = uicontrol(app.Figure,...          % Symbol edit box
                'Style','edit','Position',[75 20 50 20 ],...
                'String',app.TickerSymbol,...
                'Callback', @app.symbolUpdateCallback);

            prices = NaN*ones(1,app.NumValues) ;               % Initialize prices
            app.Line = plot(app.Axis,prices,'Marker',...
                '.','LineStyle','-');
            ylabel(app.Axis,'Stock value ($)') ;
            set(app.Axis,'XTickLabel','') ;
            title(app.Axis,['Stock Price: ' app.TickerSymbol])

            app.Timer = timer;                                 % Create timer
            app.Timer.ExecutionMode = 'fixedRate' ;
            app.Timer.Period = app.TimerUpdateRate ;
            app.Timer.TimerFcn = @app.valueUpdateCallback ;
            start(app.Timer) ;
        end

        function closeApp(app,hObject,eventdata)
        % This function runs when the app is closed
            try
                stop(app.Timer)
                delete(app.Timer)
            end
            delete(app.Figure)
        end

        function symbolUpdateCallback(app,hObject,eventdata)
        % This function runs when ticker changes in edit box
            set(app.Line,'YData',NaN*ones(1,app.NumValues));
            app.TickerSymbol = get(app.TickerEdit,'String');
            [price,name] = getQuote(app) ;
            if price == 0
                warndlg(['Ticker symbol ' app.TickerSymbol ' not found'])
            end

            title(app.Axis,['Stock Price: ' name])
            valueUpdateCallback(app)
        end

        function valueUpdateCallback(app,hObject,eventdata)
        % This function runs when the timer updates
            StockSymbol = app.TickerSymbol;
            try
                price = getQuote(app) ;
            catch
                errordlg(['Could not retrieve price for ' StockSymbol])
            end
            yvalues = get(app.Line,'YData');                   % Update the plot
            yvalues = [yvalues(2:end) price];
            set(app.Line,'YData',yvalues)
        end

        function [value,name] = getQuote(app)
            % getQuote(app)  Get stock quote using simple Yahoo finance API.
            % See http://www.gummy-stuff.org/Yahoo-data.htm for details
            url = sprintf('http://finance.yahoo.com/d/quotes.csv?s=%s&f=nl1',...
                app.TickerSymbol);
            s = urlread(url);
            [name,remain] = strtok(s,'"');
            value = str2num(strtok(remain,'",'));
        end

    end
end                                                      % End of class definition
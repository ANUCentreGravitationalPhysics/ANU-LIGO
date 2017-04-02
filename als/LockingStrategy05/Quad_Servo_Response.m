%% Quad Servo Feedback to both the ETMs (out of phase). This
%% requires the Quad response and its servo, for now I have a single
%% pendulum replacing the Quad...
    %    
    % Quad LSC input to Differential Mode error signal
    hdl = figure(13)
    a = SYS(19,16);
    b = SYS(17,16);
    G = mybodesys(b/a,f);
    %
    subplot(211)
    semilogx(f,20*log10(abs(G)))
    tt= title('Quad LSC input to Differentail Mode Error Signal TF', 'FontSize', 16);
    ylabel('Mag [dB]', 'FontSize', 14);
    axis([min(f) max(f) floor(min(log10(abs(G))))*20 ceil(max(log10(abs(G))))*20])
    grid on
    subplot(212)
    semilogx(f,angle(G)*180/pi)
    ylabel('Phase [deg]', 'FontSize', 14);
    xlabel('Freq [Hz]', 'FontSize', 14);
    axis([min(f) max(f) 90*floor(min(angle(G))/90) 90*ceil(max(angle(G))/90)])
%    set(gca,'YTick',[-1800:90:1800])
    grid on
    if save_figure == 1
        orient(hdl, 'landscape');
        print('-dpdf', ['sim/',get(tt,'string'), '.pdf']);
    end

    
    % Diff Mode Servo Open Loop response
    hdl= figure(14)
    a = SYS(17,15);
    b = SYS(19,15);
    G = mybodesys(b/a,f);
    %
    subplot(211)
    semilogx(f,20*log10(abs(G)))
    tt= title('Differential Mode Open Loop TF', 'FontSize',16);
    ylabel('Mag [dB]', 'FontSize', 14)
    axis([min(f) max(f) floor(min(log10(abs(G))))*20 ceil(max(log10(abs(G))))*20])
    grid on
    subplot(212)
    semilogx(f,angle(G)*180/pi)
    ylabel('Phase [deg]', 'FontSize', 14);
    xlabel('Freq [Hz]', 'FontSize', 14);
    axis([min(f) max(f) 90*floor(min(angle(G))/90) 90*ceil(max(angle(G))/90)])
%    set(gca,'YTick',[-1800:90:1800])
    grid on
    if save_figure == 1
        orient(hdl, 'landscape');
        print('-dpdf', ['sim/',get(tt,'string'), '.pdf']);
    end
    
    % Differentail Mode Close Loop Response
    hdl = figure(15)
    a = SYS(17,15);
    b = SYS(18,15);
    G = mybodesys(b/a,f);
    %
    subplot(211)
    semilogx(f,20*log10(abs(G)))
    tt= title('Differential Mode Closed Loop TF', 'FontSize', 16);
    ylabel('Mag [dB]', 'FontSize', 14);
    axis([min(f) max(f) floor(min(log10(abs(G))))*20 ceil(max(log10(abs(G))))*20])
    grid on
    subplot(212)
    semilogx(f,angle(G)*180/pi)
    ylabel('Phase [deg]', 'FontSize', 14);
    xlabel('Freq [Hz]', 'FontSize', 14);
    axis([min(f) max(f) 90*floor(min(angle(G))/90) 90*ceil(max(angle(G))/90)])
%    set(gca,'YTick',[-1800:90:1800])
    grid on
    if save_figure == 1
        orient(hdl, 'landscape');
        print('-dpdf', ['sim/',get(tt,'string'), '.pdf']);
    end
%%
    % Differential Mode Supression Response
    hdl = figure(16)
    G = mybodesys(SYS(17,15),f);
    %
    subplot(211)
    semilogx(f,20*log10(abs(G)))
    tt= title('Differential Suppression Response', 'FontSize', 16);
    ylabel('Mag [dB]', 'FontSize', 14);
    axis([min(f) max(f) floor(min(log10(abs(G))))*20 ceil(max(log10(abs(G))))*20])
    grid on
    subplot(212)
    semilogx(f,angle(G)*180/pi)
    ylabel('Phase [deg]', 'FontSize', 14);
    xlabel('Freq [Hz]', 'FontSize', 14);
    axis([min(f) max(f) 90*floor(min(angle(G))/90) 90*ceil(max(angle(G))/90)])
%    set(gca,'YTick',[-1800:90:1800])
    grid on
    if save_figure == 1
        orient(hdl, 'landscape');
        print('-dpdf', ['sim/',get(tt,'string'), '.pdf']);
    end
    

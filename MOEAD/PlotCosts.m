function PlotCosts(EP)

    EPC = [EP.Cost];
    plot(EPC(1, :), EPC(2, :), 'x');
    xlabel('1^{st} Objective');
    ylabel('2^{nd} Objective');
    grid on;

end
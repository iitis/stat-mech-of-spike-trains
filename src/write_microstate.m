function write_microstate(ms, fname)
    % ---- wiersze=neurony, kolumny=czas ----
    writematrix(ms, fname, 'Delimiter', ' ');
end

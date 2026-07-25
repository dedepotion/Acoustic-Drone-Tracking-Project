% Filtre passe-bande zero-phase applique aux 4 pistes du reseau de micros
% avant le calcul des tau_ij (xcorr) dans user.mlx / cercles.mlx
%
% Objectif : retirer le bruit de vent basse frequence (0-150 Hz, tres
% energetique) et le bruit haute frequence (>5 kHz) qui masquent la
% signature acoustique du drone (helices/moteurs), sans introduire de
% dephasage entre pistes (crucial pour les mesures de delai inter-micros).

lowcut  = 500;    % Hz - adjust to fit the drone's signature
highcut = 5000;   % Hz
ordre   = 4;

for i = 1:4
    infile  = sprintf('Audio Track-%d.wav', i);
    outfile = sprintf('Audiofiltre%d.wav', i);

    [x, fs] = audioread(infile);
    x = x - mean(x);

    nyq = fs/2;
    [b, a] = butter(ordre, [lowcut highcut]/nyq, 'bandpass');

    % filtfilt = filtrage a phase nulle : ne decale pas le signal dans le
    % temps, donc ne fausse pas les tau calcules par xcorr ensuite.
    y = filtfilt(b, a, x);

    % normalisation pour eviter l'ecretage, avec un peu de marge
    y = y / max(abs(y)) * 0.95;

    audiowrite(outfile, y, fs);
    fprintf('%s -> %s (fs=%d Hz)\n', infile, outfile, fs);
end

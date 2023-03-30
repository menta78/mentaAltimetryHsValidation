clear all
close all

workdir = 'data\tidalModelPairsNew\';

nstations = 612;

Model = [];
Obs = [];
Time = [];

%adat_scatter(1:nstations) = struct('stname', {}, 'stid', {}, 'lon', {}, 'lat', {}, 'time', {}, 'model', {}, 'obs', {});
adat_scatter = struct('stname', {}, 'stid', {}, 'lon', {}, 'lat', {}, 'time', {}, 'model', {}, 'obs', {});

for station=0:10
    for year=1985:2022
        fl = sprintf('\\%d\\TidalGauge_GESLA_station_%d_%d0101_%d0101.npy', year, station, year, year+1);
        pthfl = fullfile(workdir, fl);
        if ~exist(pthfl, 'file')
            continue;
        end
        %pthFiles = cat(1, pthFiles, pthfl);
        satdts = readNPY(pthfl);
        sshsat = satdts(:, 1);
        sshmdl = satdts(:, 2);
        cnd = (sshsat > -100) & (sshmdl > -100);
        satdts = satdts(cnd, :);
        cnd = (satdts(:, 1) < 100) & (satdts(:, 2) < 100);
        satdts = satdts(cnd, :);

        if size(satdts, 1) > 0
            model = satdts(:, 1);
            obs = satdts(:, 2);
            lon = satdts(1, 3);
            lat = satdts(1, 4);
            time = satdts(:, 7);
        end
        
        Time = vertcat(Time, time);
        Model = vertcat(Model, model);
        Obs = vertcat(Obs, obs);
    end
    adat_scatter(station+1).stid = station + 1;
    adat_scatter(station+1).lon = lon(1);
    adat_scatter(station+1).lat = lat(1);
    adat_scatter(station+1).time = Time;
    adat_scatter(station+1).model = Model;
    adat_scatter(station+1).obs = Obs;
end

%save('adat_tidalGaugesPairs.mat', 'adat_scatter');
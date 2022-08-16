% some folders :
%'/home/praz/Documents/MASC/sample_snow_20150620/sample_cool/processed_05-Feb-2016/DATA/GOOD';
%'/media/praz/Masc-Data/APRES3_2015/PROCESSED_20160119/DATA/GOOD';


function [X,Xlab,Xname,Xt] = load_processed_data(dirname,t_str_start,t_str_stop,feat_vec)

    fprintf('Load processed data...');

    tmin = datenum(t_str_start,'yyyymmddHHMMSS');
    tmax = datenum(t_str_stop,'yyyymmddHHMMSS');

    file_list = dir(fullfile(dirname,'20*.mat'));
    % for Massimo data
    if isempty(file_list)
        file_list = dir(fullfile(dirname,'ICE*.mat'));
    end
    % for Vanderbilt triplets
    if isempty(file_list)
        file_list = dir(fullfile(dirname,'FLAKE*.mat'));
    end
    file_only_list = {file_list.name};
    file_list = fullfile(dirname,file_only_list);

    % loop over the snowflakes and save each feature in a vector of the main
    % structure data
    n_flakes = 0;
    X = [];
    Xlab = {};
    Xname = {};
    Xt = [];
    k = 0;
    for i=1:length(file_list)

        load(file_list{i});
        % load only pictures within the time interval
        if roi.tnum >= tmin && roi.tnum <= tmax

            k = k + 1;
            j = 1;
            n_flakes = n_flakes + 1;
            
            Xname{k} = file_only_list{i};
            %disp(Xname{i});
            Xt(k) = roi.tnum;

            % dimension related features 
            X(k,j) = roi.area; Xlab{j} = 'area'; j=j+1; %1
            %X(i,j) = roi.Dmean; Xlab{j} = 'dim mean'; j=j+1;
            X(k,j) = roi.Dmax; Xlab{j} = 'Dmax'; j=j+1; %2
            X(k,j) = roi.width; Xlab{j} = 'width'; j=j+1; %3
            X(k,j) = roi.height; Xlab{j} = 'height'; j=j+1; %4
            X(k,j) = roi.perim; Xlab{j} = 'perim'; j=j+1; %5
            X(k,j) = roi.eq_radius; Xlab{j} = 'eq radius'; j=j+1; %6
            X(k,j) = roi.area_porous; Xlab{j} = 'porous area'; j=j+1; %7
            X(k,j) = (roi.area - roi.area_porous)/roi.area; Xlab{j} = 'porous ratio'; j=j+1; %8

            % ellipse related features
            X(k,j) = roi.E.a; Xlab{j} = 'ellipse fit A dim'; j=j+1; %9
            X(k,j) = roi.E.b; Xlab{j} = 'ellipse fit B dim'; j=j+1; %10
            X(k,j) = roi.E.theta; Xlab{j} = 'orientation'; j=j+1; %11
            X(k,j) = roi.E_in.a; Xlab{j} = 'ellipse in A dim'; j=j+1; %12
            X(k,j) = roi.E_in.b; Xlab{j} = 'ellipse in B dim'; j=j+1; %13
            X(k,j) = roi.E_out.a; Xlab{j} = 'ellipse out A dim'; j=j+1; %14
            X(k,j) = roi.E_out.b; Xlab{j} = 'ellipse out B dim'; j=j+1; %15
            X(k,j) = pi*roi.E_in.a*roi.E_in.b; Xlab{j} = 'ellipse in area'; j=j+1; %16
            X(k,j) = pi*roi.E_out.a*roi.E_out.b; Xlab{j} = 'ellipse out area'; j=j+1; %17
            X(k,j) = pi*roi.E.a*roi.E.b; Xlab{j} = 'ellipse fit area'; j=j+1; %18
            
            
            X(k,j) = roi.E.a/roi.E_out.a; Xlab{j} = 'ratio ellipse A fit/out'; j=j+1; %19
            X(k,j) = roi.E.b/roi.E_out.b; Xlab{j} = 'ratio ellipse B fit/out'; j=j+1; %20 
            X(k,j) = roi.E_in.a/roi.E_out.a; Xlab{j} = 'ratio ellipse A in/out'; j=j+1; %21
            X(k,j) = roi.E_in.b/roi.E_out.b; Xlab{j} = 'ratio ellipse B in/out'; j=j+1; %22
            X(k,j) = roi.E_in.a/roi.E.a; Xlab{j} = 'ratio ellipse A in/fit'; j=j+1; %23
            X(k,j) = roi.E_in.b/roi.E.b; Xlab{j} = 'ratio ellipse B in/fit'; j=j+1; %24
            X(k,j) = (roi.E_in.a*roi.E_in.b) / (roi.E_out.a*roi.E_out.b); Xlab{j} = 'ratio ellipse area in/out'; j=j+1; %25
            X(k,j) = (roi.E_in.a*roi.E_in.b) / (roi.E.a*roi.E.b); Xlab{j} = 'ratio ellipse area in/fit'; j=j+1; %26
            X(k,j) = (roi.E.a*roi.E.b) / (roi.E_out.a*roi.E_out.b); Xlab{j} = 'ratio ellipse area fit/out'; j=j+1; %27
            
            

            % intensity related features
            X(k,j) = roi.mean_intens; Xlab{j} = 'mean intens'; j=j+1; %28
            X(k,j) = roi.max_intens; Xlab{j} = 'max intens'; j=j+1; %29
            X(k,j) = roi.range_intens; Xlab{j} = 'range intens'; j=j+1; %30
            X(k,j) = roi.focus; Xlab{j} = 'focus'; j=j+1; %31
            X(k,j) = roi.area_focus; Xlab{j} = 'area focus'; j=j+1; %32
            X(k,j) = roi.area_range; Xlab{j} = 'area range'; j=j+1; %33
            X(k,j) = roi.lap; Xlab{j} = 'lap energy'; j=j+1; %34
            X(k,j) = roi.area_lap; Xlab{j} = 'area X lap'; j=j+1; %35
            X(k,j) = roi.std; Xlab{j} = 'global std'; j=j+1; %36
            X(k,j) = roi.local_std; Xlab{j} = 'avg local std 3x3'; j=j+1; %37
            X(k,j) = roi.local_std5; Xlab{j} = 'avg local std 5x5'; j=j+1; %38
            X(k,j) = roi.local_std7; Xlab{j} = 'avg local std 7x7'; j=j+1; %39
            X(k,j) = roi.contrast; Xlab{j} = 'contrast'; j=j+1; %40

            % enhanced CLAHE picture related features
            X(k,j) = roi.new.range_intens; Xlab{j} = 'bright range intens'; j=j+1; %41
            X(k,j) = roi.new.lap; Xlab{j} = 'bright lap energy'; j=j+1; %42
            X(k,j) = roi.new.area_lap; Xlab{j} = 'bright area X lap'; j=j+1; %43
            X(k,j) = roi.new.std; Xlab{j} = 'bright global std'; j=j+1; %44
            X(k,j) = roi.new.local_std; Xlab{j} = 'bright local std 3x3'; j=j+1; %45
            X(k,j) = roi.new.local_std5; Xlab{j} = 'bright local std 5x5'; j=j+1; %46
            X(k,j) = roi.new.local_std7; Xlab{j} = 'bright local std 7x7'; j=j+1; %47
            X(k,j) = roi.new.contrast; Xlab{j} = 'bright contrast'; j=j+1; %48

            % Skeleton
            X(k,j) = roi.skel.p_ratio; Xlab{j} = 'skel/perim ratio'; j=j+1; %49
            X(k,j) = roi.skel.A_ratio; Xlab{j} = 'skel/area ratio'; j=j+1; %50
            if ~isnan(roi.skel.N_ends)
                X(k,j) = roi.skel.N_ends; Xlab{j} = 'skel N ends'; j=j+1; %51  
            else
                X(k,j) = 0; Xlab{j} = 'skel N ends'; j=j+1; %51  
            end
            X(k,j) = roi.skel.N_junctions; Xlab{j} = 'skel N junctions'; j=j+1; %52

            % Haralick 49-52
            X(k,j) = roi.H.Contrast; Xlab{j} = 'Haralick contrast'; j=j+1; %53
            X(k,j) = roi.H.Correlation; Xlab{j} = 'Haralick correlation'; j=j+1; %54
            X(k,j) = roi.H.Energy; Xlab{j} = 'Haralick energy'; j=j+1; %55
            X(k,j) = roi.H.Homogeneity; Xlab{j} = 'Haralick homogeneity'; j=j+1; %56

            % Fractal
            X(k,j) = roi.F; Xlab{j} = 'boxcounting fractal dim'; j=j+1; %57
            X(k,j) = roi.F_jac; Xlab{j} = 'theoretical fractal dim'; j=j+1; %58

            % Hull
            X(k,j) = roi.hull.solidity; Xlab{j} = 'solidity'; j=j+1; %59
            X(k,j) = roi.hull.convexity; Xlab{j} = 'convexity'; j=j+1; %60
            X(k,j) = length(roi.hull.xh); Xlab{j} = 'Hull nb angles'; j=j+1; %61
            
            
            % Others
            X(k,j) = max(1,roi.complex); Xlab{j} = 'complexity'; j=j+1; %62
            X(k,j) = roi.roundness; Xlab{j} = 'roundness'; j=j+1; %63
            X(k,j) = roi.compactness; Xlab{j} = 'compactness'; j=j+1; %64
            if roi.E.a >= roi.E.b
                X(k,j) = roi.E.b/roi.E.a; Xlab{j} = 'aspec ratio'; j=j+1; %65
                X(k,j) = sqrt(1-roi.E.b/roi.E.a); Xlab{j} = 'eccentricity'; j=j+1; %66
            else
                X(k,j) = 1; Xlab{j} = 'aspect ratio'; j=j+1; %65
                X(k,j) = 0; Xlab{j} = 'eccentricity'; j=j+1; %66
            end
            X(k,j) = roi.xhi; Xlab{j} = 'xi'; j=j+1; %67
            
            if roi.nb_holes > 0
                X(k,j) = 1; Xlab{j} = 'has holes'; j=j+1; %68
            else
                X(k,j) = 0; Xlab{j} = 'has holes'; j=j+1; %68
            end
            X(k,j) = roi.nb_holes; Xlab{j} = 'N holes'; j=j+1; %69
            
            % New descriptors based on min rectangle box
            X(k,j) = roi.Rect.A_ratio; Xlab{j} = 'rectangularity'; j=j+1; %70
            X(k,j) = roi.Rect.p_ratio; Xlab{j} = 'rect perim ratio'; j=j+1; %71
            X(k,j) = roi.Rect.aspect_ratio; Xlab{j} = 'rect aspect ratio'; j=j+1; %72
            X(k,j) = roi.Rect.eccentricity; Xlab{j} = 'rect eccentricity'; j=j+1; %73
            
            % New descriptor based on circumscribed circle
            X(k,j) = (2*pi*roi.C_out.r)/roi.perim; Xlab{j} = 'circle out perim ratio'; j=j+1; %74
            
            
%             if ~isempty(roi.fallspeed)
%                 X(i,j) = roi.fallspeed; Xlab{j} = 'fallspeed'; j=j+1;
%             else
%                 X(i,j) = NaN; Xlab{j} = 'fallspeed'; j=j+1;
%             end
            
            % 3 new descriptors never used so far
            X(k,j) = roi.wavs; Xlab{j} = 'wavelet desc'; j=j+1; %75
            X(k,j) = roi.hist_entropy; Xlab{j} = 'hist_entropy'; j=j+1; %76
            X(k,j) = max(1,roi.complex) * (1 + roi.range_intens); Xlab{j} = 'chi Garrett wrong'; j=j+1; %77
            X(k,j) = max(1,roi.complex) * (1 + roi.local_std7); Xlab{j} = 'chi Garrett true'; j=j+1; %78
            X(k,j) = roi.local_std; Xlab{j} = 'std_intens'; j=j+1; %79
            X(k,j) = roi.complex; Xlab{j} = 'complex'; j=j+1; %80
            
            if 1
            
            % New descriptors which are probably useless
            X(k,j) = roi.Sym.P0; Xlab{j} = 'fft P0'; j=j+1; %81
            X(k,j) = roi.Sym.P1; Xlab{j} = 'fft P1'; j=j+1; %82
            X(k,j) = roi.Sym.P2; Xlab{j} = 'fft P2'; j=j+1; %83
            X(k,j) = roi.Sym.P3; Xlab{j} = 'fft P3'; j=j+1; %84
            X(k,j) = roi.Sym.P4; Xlab{j} = 'fft P4'; j=j+1; %85
            X(k,j) = roi.Sym.P5; Xlab{j} = 'fft P5'; j=j+1; %86
            X(k,j) = roi.Sym.P6; Xlab{j} = 'fft P6'; j=j+1; %87
            X(k,j) = roi.Sym.P7; Xlab{j} = 'fft P7'; j=j+1; %88
            X(k,j) = roi.Sym.P8; Xlab{j} = 'fft P8'; j=j+1; %89
            X(k,j) = roi.Sym.P9; Xlab{j} = 'fft P9'; j=j+1; %90
            X(k,j) = roi.Sym.P10; Xlab{j} = 'fft P10'; j=j+1; %91
            X(k,j) = roi.Sym.mean; Xlab{j} = 'Sym mean'; j=j+1; %92
            X(k,j) = roi.Sym.std; Xlab{j} = 'Sym std'; j=j+1; %93
            X(k,j) = roi.Sym.std/roi.Sym.mean; Xlab{j} = 'Sym std/mean'; j=j+1; %94
            X(k,j) = roi.Sym.P6/max([roi.Sym.P1,roi.Sym.P2,roi.Sym.P3,roi.Sym.P4,roi.Sym.P5,roi.Sym.P6,roi.Sym.P7,roi.Sym.P8,roi.Sym.P9,roi.Sym.P10]); Xlab{j} = 'fft P6/Pmax'; j=j+1; %95
            [~,idx_max] = max([roi.Sym.P1,roi.Sym.P2,roi.Sym.P3,roi.Sym.P4,roi.Sym.P5,roi.Sym.P6,roi.Sym.P7,roi.Sym.P8,roi.Sym.P9,roi.Sym.P10]); 
            X(k,j) = idx_max-1; Xlab{j} = 'fft P# max'; j=j+1; %96
 
            % New descriptors based on focal measure
%             X(i,j) = roi.fmeas.acmo; Xlab{j} = 'ACMO'; j=j+1; %97
%             X(i,j) = roi.fmeas.bren; Xlab{j} = 'BREN'; j=j+1; %98
%             X(i,j) = roi.fmeas.cont; Xlab{j} = 'CONT'; j=j+1; %99
%             X(i,j) = roi.fmeas.curv; Xlab{j} = 'CURV'; j=j+1; %100
%             X(i,j) = roi.fmeas.dctr; Xlab{j} = 'DCTR'; j=j+1; %101 - NANANANAN
%             X(i,j) = roi.fmeas.gder; Xlab{j} = 'GDER'; j=j+1; %102
%             X(i,j) = roi.fmeas.glva; Xlab{j} = 'GLVA'; j=j+1; %103
%             X(i,j) = roi.fmeas.gllv; Xlab{j} = 'GLLV'; j=j+1; %104
%             X(i,j) = roi.fmeas.glvn; Xlab{j} = 'GLVN'; j=j+1; %105
%             X(i,j) = roi.fmeas.grae; Xlab{j} = 'GRAE'; j=j+1; %106
%             X(i,j) = roi.fmeas.grat; Xlab{j} = 'GRAT'; j=j+1; %107
%             X(i,j) = roi.fmeas.gras; Xlab{j} = 'GRAS'; j=j+1; %108
%             X(i,j) = roi.fmeas.helm; Xlab{j} = 'HELM'; j=j+1; %109
%             X(i,j) = roi.fmeas.hisr; Xlab{j} = 'HISR'; j=j+1; %110
%             X(i,j) = roi.fmeas.lape; Xlab{j} = 'LAPE'; j=j+1; %111
%             X(i,j) = roi.fmeas.lapv; Xlab{j} = 'LAPV'; j=j+1; %112
%             X(i,j) = roi.fmeas.lapd; Xlab{j} = 'LAPD'; j=j+1; %113
%             X(i,j) = roi.fmeas.sfil; Xlab{j} = 'SFIL'; j=j+1; %114
%             X(i,j) = roi.fmeas.sfrq; Xlab{j} = 'SFRQ'; j=j+1; %115
%             X(i,j) = roi.fmeas.teng; Xlab{j} = 'TENG'; j=j+1; %116
%             X(i,j) = roi.fmeas.tenv; Xlab{j} = 'TENV'; j=j+1; %117
%             X(i,j) = roi.fmeas.vola; Xlab{j} = 'VOLA'; j=j+1; %118
%             X(i,j) = roi.fmeas.wavv; Xlab{j} = 'WAVV'; j=j+1; %119
%             X(i,j) = roi.fmeas.wavr; Xlab{j} = 'WAVR'; j=j+1; %120
            
            
            % 77 here
            
            % Full set of Haralick features (20)
            %X(i,j) = mean(roi.H_new.autoc); Xlab{j} = 'Haralick autocorrelation'; j=j+1;
            %X(i,j) = mean(roi.H_new.contr); Xlab{j} = 'Haralick contrast'; j=j+1;
            %X(i,j) = mean(roi.H_new.corrm); Xlab{j} = 'Haralick correlation'; j=j+1;
            %X(i,j) = mean(roi.H_new.cprom); Xlab{j} = 'Haralick cluster prominence'; j=j+1;
%             X(i,j) = mean(roi.H_new.cshad); Xlab{j} = 'Haralick cluster shade'; j=j+1;
%             X(i,j) = mean(roi.H_new.dissi); Xlab{j} = 'Haralick dissimilarity'; j=j+1;
%             X(i,j) = mean(roi.H_new.energ); Xlab{j} = 'Haralick energy'; j=j+1;
%             X(i,j) = mean(roi.H_new.entro); Xlab{j} = 'Haralick entropy'; j=j+1;
%             X(i,j) = mean(roi.H_new.homom); Xlab{j} = 'Haralick homogeneity'; j=j+1;
%             X(i,j) = mean(roi.H_new.maxpr); Xlab{j} = 'Haralick max prob'; j=j+1;
%             X(i,j) = mean(roi.H_new.sosvh); Xlab{j} = 'Haralick sum of squares variance'; j=j+1;
%             X(i,j) = mean(roi.H_new.savgh); Xlab{j} = 'Haralick sum average'; j=j+1;
%             X(i,j) = mean(roi.H_new.svarh); Xlab{j} = 'Haralick sum variance'; j=j+1;
%             X(i,j) = mean(roi.H_new.senth); Xlab{j} = 'Haralick sum entropy'; j=j+1;
%             X(i,j) = mean(roi.H_new.dvarh); Xlab{j} = 'Haralick difference variance'; j=j+1;
%             X(i,j) = mean(roi.H_new.denth); Xlab{j} = 'Haralick difference entropy'; j=j+1;
%             X(i,j) = mean(roi.H_new.inf1h); Xlab{j} = 'Haralick info measure of corr 1'; j=j+1;
%             X(i,j) = mean(roi.H_new.inf2h); Xlab{j} = 'Haralick info measure of corr 2'; j=j+1;
%             X(i,j) = mean(roi.H_new.indnc); Xlab{j} = 'Haralick inverse diff normalized (INN)'; j=j+1;
%             X(i,j) = mean(roi.H_new.idmnc); Xlab{j} = 'Haralick inverse diff moment normalized'; j=j+1;
            
            % descriptors based on new skeleton
%             try
%                 X(i,j) = roi.skel2.max_dist; Xlab{j} = 'skel2 max dist'; j=j+1;
%                 X(i,j) = roi.skel2.min_dist; Xlab{j} = 'skel2 min dist'; j=j+1;
%                 X(i,j) = roi.skel2.mean_dist; Xlab{j} = 'skel2 mean dist'; j=j+1;
%                 X(i,j) = roi.skel2.med_dist; Xlab{j} = 'skel2 med dist'; j=j+1;
%                 X(i,j) = roi.skel2.std_dist; Xlab{j} = 'skel2 std dist'; j=j+1;
% 
%                 X(i,j) = roi.skel2.max_interangle; Xlab{j} = 'skel2 max interangle'; j=j+1;
%                 X(i,j) = roi.skel2.min_interangle; Xlab{j} = 'skel2 min interangle'; j=j+1;
%                 X(i,j) = roi.skel2.med_interangle; Xlab{j} = 'skel2 med interangle'; j=j+1;
%                 X(i,j) = roi.skel2.mean_interangle; Xlab{j} = 'skel2 mean interangle'; j=j+1;
%                 X(i,j) = roi.skel2.std_interangle; Xlab{j} = 'skel2 std interangle'; j=j+1;
%             
%                 X(i,j) = roi.skel2.F; Xlab{j} = 'skel2 fractal dim'; j=j+1;
%                 
%             catch
%                 X(i,j) = 0; Xlab{j} = 'skel2 max dist'; j=j+1;
%                 X(i,j) = 0; Xlab{j} = 'skel2 min dist'; j=j+1;
%                 X(i,j) = 0; Xlab{j} = 'skel2 mean dist'; j=j+1;
%                 X(i,j) = 0; Xlab{j} = 'skel2 med dist'; j=j+1;
%                 X(i,j) = 0; Xlab{j} = 'skel2 std dist'; j=j+1;
% 
%                 X(i,j) = 0; Xlab{j} = 'skel2 max interangle'; j=j+1;
%                 X(i,j) = 0; Xlab{j} = 'skel2 min interangle'; j=j+1;
%                 X(i,j) = 0; Xlab{j} = 'skel2 med interangle'; j=j+1;
%                 X(i,j) = 0; Xlab{j} = 'skel2 mean interangle'; j=j+1;
%                 X(i,j) = 0; Xlab{j} = 'skel2 std interangle'; j=j+1;
%             
%                 X(i,j) = 0; Xlab{j} = 'skel2 fractal dim'; j=j+1;
%                 fprintf('Warning: in ROI %s:\n',file_only_list{i});
%                 fprintf('skel2 attributes have been set to 0 because they were not defined in the ROI.\n');
%             end
                
            %X(i,j) = roi.skel2.F_4box; Xlabl{j} = 'skel2 fractal dim 4 boxes'; j=j+1;
            
            %X(i,j) = roi.skel2.p_ratio; Xlab{j} = 'skel2 p_ratio'; j=j+1;
            %X(i,j) = roi.skel2.A_ratio; Xlab{j} = 'skel2 A_ratio'; j=j+1;
           
            %X(i,j) = roi.skel2.N_ends; Xlab{j} = 'skel2 N ends'; j=j+1;
            %X(i,j) = roi.skel2.N_junctions; Xlab{j} = 'skel2 N junctions'; j=j+1;
            
            % More shit
            %disp(roi.glcm_vec32');
            %disp(length(roi.glcm_vec32));
            %X(i,j:j+length(roi.glcm_vec32)-1) = roi.glcm_vec32'; %j=j+length(roi.glcm_vec32);
            
            % New shit
            %X(i,j) = roi.Hu.h1; Xlab{j} = 'Hu #1'; j=j+1;
            %X(i,j) = roi.Hu.h2; Xlab{j} = 'Hu #2'; j=j+1;
            %X(i,j) = roi.Hu.h3; Xlab{j} = 'Hu #3'; j=j+1;
            %X(i,j) = roi.Hu.h4; Xlab{j} = 'Hu #4'; j=j+1;
            %X(i,j) = roi.Hu.h5; Xlab{j} = 'Hu #5'; j=j+1;
            %X(i,j) = roi.Hu.h6; Xlab{j} = 'Hu #6'; j=j+1;
            %X(i,j) = roi.Hu.h7; Xlab{j} = 'Hu #7'; j=j+1;
            %X(i,j) = roi.Hu.h8; Xlab{j} = 'Hu #8'; j=j+1;
            
           

            end

        end

    end
    
    % 4th argument is list of features wanted
    if nargin==4 && ~isempty(X)

        X = X(:,feat_vec);
        Xlab = Xlab(feat_vec);

    end
    
    
    fprintf('   Done! %u images found and processed. %u descriptors computed. \n',n_flakes,size(X,2));
    Xlab = Xlab';
    Xname = Xname';
    Xt = Xt';
    
end

% === Core code to process MASC data =====================================
%
% Load one image and process it :
%
% 1) masking (see function help)
% 2) creating a B&W mask
% 3) find the "best" ROI (see function roi_detection.m)
% 4) crop around it, clean leftover pixels
% 5) compute some features (see the code)
% 6) plot some illustrations if desired
% 7) save the cropped image in a .png file and a structure of features in a
% .mat structure if desired
%
% Adapted from Tim Garrett, University of Utah. 
% Author : Christophe Praz (christophe.praz@epfl.ch)
% Last update : October 2015
% ========================================================================
function [tt,flag] = MASC_picture_process(j,current_dir_list,pic_list,label,process)

            % Load img if possible
            t_loading = tic;
            try

                flake_filepath = fullfile(current_dir_list,pic_list.files{j});
                flake_filename = pic_list.files{j};
                %flake_imfinfo = imfinfo(flake_filepath);
                flake_data = imread(flake_filepath);  
                if length(size(flake_data)) > 2
                    flake_data = rgb2gray(flake_data);
                end
                flake_id = pic_list.id(j);
                flake_fallspeed = pic_list.fallspeed(pic_list.fallid == flake_id);
                flake_cam = pic_list.cam(j);

                % resolution [um]
                %flake_XRes = cam.fovmat(flake_cam+1)/flake_imfinfo.Width*1000;
                %flake_YRes = cam.fovmat(flake_cam+1)/flake_imfinfo.Height*1000;


            catch

                fprintf('Warning : image %s is skipped (couldn''t load the image)\n',flake_filename);
                flag = -1;
                tt = [];
                return;

            end   
            tt.loading = toc(t_loading);

            % Remove clutter
            t_clutter = tic;
            data = masking(flake_data,flake_cam,process);
            tt.clutter = toc(t_clutter);
            
            % Simply create a B&W mask 
            t_edging = tic;
            data_eroded = false(size(data));
            data_eroded(data>0) = 1;
            data_eroded = imfill(data_eroded,'holes');
            tt.edging = toc(t_edging);
            
            % Detecting the ROIs
            t_roiying = tic; 
            [all_roi,idx_best,area_focus_ratio,flag,status] = roi_detection(data,data_eroded,process,flake_cam,true);
            %roi = roi_detection(data,data_eroded_new,process,flake_cam);               
            tt.roiying = toc(t_roiying);
            
            % If no ROI found at all
            if isempty(all_roi)

                %fprintf('WARNING : no ROI candidate found in %s \n', flake_filename);
                flag = 0;
                tt.feature = 0;
                tt.plotting = 0;
                tt.saving = 0;
                return;  
            end
                
            % this is just for Vanderbilt triplet data (already cropped)
            
%             elseif numel(all_roi) > 1
%                 
%                 fprintf('WARNING : %u ROIs found in %s \n',numel(all_roi),flake_filename);
% 
%             end
            
            % Otherwise we continue and compute features on the best ROI
            t_feature = tic;
            
            roi = process_basic_descriptors(data,all_roi(idx_best),process);
            % additional descriptors not based on the image itself
            roi.area_focus_ratio = area_focus_ratio;
            roi.flag = flag;
            roi.status = status;
            roi.n_roi = length(all_roi);
            roi.name = flake_filename;
            roi.id = flake_id;
            roi.cam = flake_cam;
            roi.tnum = pic_list.time_num(j);
            roi.fallspeed = flake_fallspeed;
            
            tt.feature = toc(t_feature);
           

            % the flag is used in MASC_process to count good/bad images
            if strcmp(roi.flag,'GOOD')
                roi.status = 'good detection.';
                flag = 2;
            else 
                flag = 0;
            end
            
            % Plot generation (illustration if desired)
            t_plotting = tic;
            if process.generate_figs

%                 fig1=figure('Visible',process.display_figs,'Renderer','painters');%,'Position',[100 100 700 1300]);     
%                 subplot(5,4,1:12);
%                 imshow(flake_data);
%                 %title(sprintf('v_{f} = %1.2f [m/s]     T = %1.1f [C]     RH = %1.2f [%%]',roi.fallspeed,roi.thygan.T,roi.thygan.RH));
%                 subplot(5,4,[13 17]);
%                 imshow(roi.data);
%                 subplot(5,4,[14 18]);
%                 imshow(roi.new.data);
%                 h_text = subplot(5,4,[16 20]);
%                 xl = xlim(h_text); 
%                 xPos = xl(1) + diff(xl) / 2; 
%                 yl = ylim(h_text); 
%                 yPos = yl(1) + diff(yl) / 2; 
%                 plot_text = text(xPos, yPos, sprintf('dim : %2.0f \n bright. : %1.2f \n range : %1.3f \n new range : %1.3f \n contrast : %1.2f \n new contrast : %1.2f \n std : %1.2f \n new std : %1.2f \n local std : %1.2f \n new local std : %1.2f \n lap : %1.2f \n new lap : %1.2f \n complex : %1.2f',...
%                     double(roi.width+roi.height)/2,roi.mean_intens,roi.range_intens,roi.new.range_intens,roi.contrast,roi.new.contrast,roi.std,roi.new.std,roi.local_std,roi.new.local_std,roi.lap,roi.new.lap,roi.complex), 'Parent', h_text);
%                 set(plot_text, 'HorizontalAlignment', 'center','FontSize',12,'FontWeight','demi');
%                 set(h_text,'visible','off');           

                fig2=figure('Visible',process.display_figs);
                colormap(gray(255));
                hold on; axis equal; box on;
                image(double(roi.data));
                set(gca,'YDir','reverse');  
                
                t = linspace(0,2*pi);
                theta = -roi.E.theta*pi/180;
                xt1 = roi.E.X0 + cos(theta).*roi.E.a.*cos(t) - sin(theta).*roi.E.b.*sin(t);
                yt1 = roi.E.Y0 + sin(theta).*roi.E.a.*cos(t) + cos(theta).*roi.E.b.*sin(t);
                xt2 = roi.E_out.X0 + cos(theta).*roi.E_out.a.*cos(t) - sin(theta).*roi.E_out.b.*sin(t);
                yt2 = roi.E_out.Y0 + sin(theta).*roi.E_out.a.*cos(t) + cos(theta).*roi.E_out.b.*sin(t);
                xt3 = roi.E_in.X0 + cos(theta).*roi.E_in.a.*cos(t) - sin(theta).*roi.E_in.b.*sin(t);
                yt3 = roi.E_in.Y0 + sin(theta).*roi.E_in.a.*cos(t) + cos(theta).*roi.E_in.b.*sin(t);
                xt4 = roi.C_out.X0 + roi.C_out.r .*cos(t);
                yt4 = roi.C_out.Y0 + roi.C_out.r .*sin(t);

                plot(roi.x_perim,roi.y_perim,'y-');
                plot(xt1,yt1,'r-','linewidth',2);
                plot(xt2,yt2,'c-','linewidth',2);
                plot(roi.hull.xh,roi.hull.yh,'c--');
                plot(xt3,yt3,'g-','linewidth',2);
                plot(xt4,yt4,'b-','linewidth',2);
                plot(roi.E.X0,roi.E.Y0,'rx');
                plot(roi.E_out.X0,roi.E_out.Y0,'co');
                plot(roi.E_in.X0,roi.E_in.Y0,'gv');
                %plot([DmaxA(1) DmaxB(1)],[DmaxA(2) DmaxB(2)],'r--');
                line(roi.Rect.rectx,roi.Rect.recty);
                xlabel('x axis [pixels]');
                ylabel('y axis [pixels]');
                title(roi.E.theta);
                set(gca,'Color','k');
                
                % - insert more figures here -
                
                drawnow();

            end
            tt.plotting = toc(t_plotting);
    
           
            % save flakes processed in label_params.outdir if desired
            t_saving = tic;

            if process.saveresults   
                
                % creation of the metadata file (first time)
                if ~exist(label.outdir,'dir')    
                    mkdir(label.outdir);
                    create_proc_params_file(label.outdir,label,process);
                    fprintf('creation of a metadata file proc_params.txt...\n');
                end
                
                date_folder = datestr(roi.tnum,'yyyy.mm.dd');
                hour_folder = datestr(roi.tnum,'HH');  
                  
                if strcmp(roi.flag,'GOOD')
                    
                    path2save = fullfile(label.outdir,date_folder,hour_folder);
                %path2save_im = fullfile(path2save,'IMAGES');
                %path2save_data = fullfile(path2save,'DATA');
                %path2save_fig = fullfile(path2save,'FIGURES');

                    if ~exist(path2save,'dir')

                        mkdir(path2save);
                        
                    end
                    %mkdir(path2save_im);
                    %mkdir(path2save_data);
                    %mkdir(path2save_fig);

                    %mkdir(fullfile(path2save_im,'GOOD'));
                    %mkdir(fullfile(path2save_im,'BAD'));
                    %mkdir(fullfile(path2save_data,'GOOD'));
                    %mkdir(fullfile(path2save_data,'BAD'));
                    %mkdir(fullfile(path2save_fig,'GOOD'));
                    %mkdir(fullfile(path2save_fig,'BAD'));
                    
                    % check if there exist a metadata file and create one
                    % if it not the case
                    %if ~exist(fullfile(path2save,'proc_params.txt'),'file')
                     %   create_proc_params_file(path2save,label,process);
                     %   fprintf('creation of a metadata file proc_params.txt...\n');
                    %end
                    
                    imwrite(roi.data,fullfile(path2save,flake_filename),'png','BitDepth', 8);
                    save_mat(fullfile(path2save),strcat(flake_filename(1:end-4),'.mat'),roi);

                elseif ~strcmp(roi.status,'no ROI detected')
                  
                    path2save = fullfile(label.outdir,'BAD',date_folder,hour_folder);
                    
                    if ~exist(path2save,'dir');

                        mkdir(path2save);
                        
                    end
                    
                    imwrite(roi.data,fullfile(path2save,flake_filename),'png','BitDepth', 8);
                    save_mat(fullfile(path2save),strcat(flake_filename(1:end-4),'.mat'),roi);
                    
                end
                
                

                % saving the image
                %fprintf('Saving %s...\n',strcat(flake_filename));
                %if strcmp(roi.flag,'GOOD')
                %    imwrite(roi.data,fullfile(path2save_im,'GOOD',flake_filename),'png','BitDepth', 8);
                %    save_mat(fullfile(path2save_data,'GOOD'),strcat(flake_filename(1:end-4),'.mat'),roi);
                %    if process.save_figs
                %        saveas(fig1,fullfile(path2save_fig,'GOOD',strcat('fig1_',flake_filename)));
                %        %print(fig1,fullfile(path2save_fig,'GOOD',strcat('fig1_',flake_filename)),'-dpng','-r80');
                %        saveas(fig2,fullfile(path2save_fig,'GOOD',strcat('fig2_',flake_filename)));
                %    end

                %elseif ~strcmp(roi.status,'no ROI detected')
                %    imwrite(roi.data,fullfile(path2save_im,'BAD',flake_filename),'png','BitDepth', 8);
                %    save_mat(fullfile(path2save_data,'BAD'),strcat(flake_filename(1:end-4),'.mat'),roi);
                %    if process.save_figs
                %        saveas(fig1,fullfile(path2save_fig,'BAD',strcat('fig1_',flake_filename)));
                %        %print(fig1,fullfile(path2save_fig,'BAD',strcat('fig1_',flake_filename)),'-dpng','-r80');
                %        saveas(fig2,fullfile(path2save_fig,'BAD',strcat('fig2_',flake_filename)));
                %    end
                %end    

            end
            
            tt.saving = toc(t_saving);
            
            % manually clear some variables because matlab is not good to
            % do it by itself when parallel processing is enabled
            if strcmp(process.display_figs,'off')
                close all;
            end
            clearvars -except tt flag;
            %clear fig1 fig2 h_text plot_text roi data data_eroded data_outlined flake_data flake_imfinfo
            % clearvars -except tt_uploading cam current_dir_list dir_list i idx_oi label N_discarded N_flake n_id N_processed pic_list process t_startprogram t_uploading tt_clutter tt_edging tt_feature tt_loading tt_plotting tt_roiying tt_saving
         
end
















% some obsolete pieces of code
%                 elseif strcmp(roi.flag,'BLURRY')
%                     imwrite(roi.data,fullfile(path2save_im,'BLURRY',flake_filename),'png','BitDepth', 8);
%                     save_mat(fullfile(path2save_data,'BLURRY'),strcat(flake_filename(1:end-4),'.mat'),roi);
%                     %saveas(fig1,fullfile(path2save_fig,'BLURRY',strcat('fig1_',flake_filename)));
%                     %print(fig1,fullfile(path2save_fig,'BLURRY',strcat('fig1_',flake_filename)),'-dpng','-r80');  
%                     %path = fullfile(path2save_fig,'BLURRY',strcat('fig1_',flake_filename));
%                     %export_fig(path,'-r65','-opengl');% -opengl -nocrop -r100











function varargout = classification_label_program(varargin)
% classification_label_program MATLAB code for classification_label_program.fig
%      classification_label_program, by itself, creates a new classification_label_program or raises the existing
%      singleton*.
%
%      H = riming_classification_label_program returns the handle to a new classification_label_program or the handle to
%      the existing singleton*.
%
%      classification_label_program('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in classification_label_program.M with the given input arguments.
%
%      classification_label_program('Property','Value',...) creates a new classification_label_program or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before classification_label_program_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to classification_label_program_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help classification_label_program

% Last Modified by GUIDE v2.5 07-Jul-2016 16:51:50

% Begin initialization code - DO NOT EDIT mfilename
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @classification_label_program_OpeningFcn, ...
                   'gui_OutputFcn',  @classification_label_program_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before classification_label_program is made visible.
function classification_label_program_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to classification_label_program (see VARARGIN)

console_Init(hObject,handles);
handles = guidata(hObject);
set(handles.figure1, 'Name', 'Snowflakes labelling');

% load a default folder 
load_default = false;
if load_default
    handles.dir_data = '/home/praz/Desktop/lte_on_enac1files/commun1/Christophe/MASC_labelling/mixed_sample_part2';
    handles.file_list = dir(fullfile(handles.dir_data,'20*.mat'));
    % for Massimo data
    if isempty(handles.file_list)
        handles.file_list = dir(fullfile(handles.dir_data,'ICE*.mat'));
    end
    handles.file_list = {handles.file_list.name}';
    handles.idx_max = length(handles.file_list);
    % randomization
    handles.idx_rnd = randperm(handles.idx_max);
    handles.file_list = handles.file_list(handles.idx_rnd);
    handles.current_idx = 1;
    set(handles.is_melting,'Value',handles.isMelting(handles.current_idx));
    load(fullfile(handles.dir_data,handles.file_list{handles.current_idx}));
    handles.current_name = roi.name;
    handles.labelId_list = zeros(length(handles.file_list),1);
    handles.labelId_list(handles.labelId_list==0) = -1;
    handles.labelName_list = cell(length(handles.file_list),1);
    handles.rimingId_list = zeros(length(handles.file_list),1);
    handles.rimingId_list(handles.rimingId_list==0) = -1;
    handles.rimingName_list = cell(length(handles.file_list),1);
    handles.isMelting = false(length(handles.file_list),1);
    %guidata(hObject,handles);
    refresh_image(handles,roi);
    load_listbox(handles);
    disp_current_folder(handles);
    total_flakes_Update(hObject,handles);
    total_labels_Update(hObject,handles);
    total_undetermined_Update(hObject,handles);
else
    handles.dir_data = pwd;
    handles.file_list = {};
    handles.labelId_list = [];
    handles.labelName_list = {};
    handles.rimingId_list = [];
    handles.rimingName_list = {};
    handles.isMelting = [];
    set(handles.is_melting,'Value',0);
    handles.current_idx = 1;
    handles.idx_max = 1;
    handles.idx_rnd = [];
    load_listbox(handles);
    console_Update(hObject,handles,'Please select a data directory by clicking on "Browse"');
    set(gca,'Xtick',[]);
    set(gca,'Ytick',[]);
    box on;
end

handles.save_labels = false;

% Choose default command line output for classification_label_program
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes classification_label_program wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = classification_label_program_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

%% class buttons

% % --- Executes on button press in c_graupel.
% function c_graupel_Callback(hObject, eventdata, handles)
% 
%     handles.labelName_list{handles.current_idx} = 'graupel';
%     handles.labelId_list(handles.current_idx) = 1;
%     total_flakes_Update(hObject,handles);
%     total_labels_Update(hObject,handles);
%     total_undetermined_Update(hObject,handles);
%     go_next_Callback(hObject, eventdata, handles);
%     
% % --- Executes on button press in c_agg.
% function c_agg_Callback(hObject, eventdata, handles)
% 
%     handles.labelName_list{handles.current_idx} = 'aggregate';
%     handles.labelId_list(handles.current_idx) = 2;
%     total_flakes_Update(hObject,handles);
%     total_labels_Update(hObject,handles);
%     total_undetermined_Update(hObject,handles);
%     go_next_Callback(hObject, eventdata, handles)
% 
% 
% % --- Executes on button press in c_melt.
% function c_melt_Callback(hObject, eventdata, handles)
% 
%     handles.labelName_list{handles.current_idx} = 'melting';
%     handles.labelId_list(handles.current_idx) = 3;
%     total_flakes_Update(hObject,handles);
%     total_labels_Update(hObject,handles);
%     total_undetermined_Update(hObject,handles);
%     go_next_Callback(hObject, eventdata, handles)
% 
% 
% % --- Executes on button press in c_small.
% function c_small_Callback(hObject, eventdata, handles)
% 
%     handles.labelName_list{handles.current_idx} = 'small';
%     handles.labelId_list(handles.current_idx) = 4;
%     total_flakes_Update(hObject,handles);
%     total_labels_Update(hObject,handles);
%     total_undetermined_Update(hObject,handles);
%     go_next_Callback(hObject, eventdata, handles)
% 
% 
% % --- Executes on button press in c_dend.
% function c_dend_Callback(hObject, eventdata, handles)
% 
%     handles.labelName_list{handles.current_idx} = 'dendrite';
%     handles.labelId_list(handles.current_idx) = 5;
%     total_flakes_Update(hObject,handles);
%     total_labels_Update(hObject,handles);
%     total_undetermined_Update(hObject,handles);
%     go_next_Callback(hObject, eventdata, handles)
% 
% 
% % --- Executes on button press in c_plate.
% function c_plate_Callback(hObject, eventdata, handles)
% 
%     handles.labelName_list{handles.current_idx} = 'plate';
%     handles.labelId_list(handles.current_idx) = 6;
%     total_flakes_Update(hObject,handles);
%     total_labels_Update(hObject,handles);
%     total_undetermined_Update(hObject,handles);
%     go_next_Callback(hObject, eventdata, handles)
% 
% 
% % --- Executes on button press in c_col.
% function c_col_Callback(hObject, eventdata, handles)
% 
%     handles.labelName_list{handles.current_idx} = 'column';
%     handles.labelId_list(handles.current_idx) = 7;
%     total_flakes_Update(hObject,handles);
%     total_labels_Update(hObject,handles);
%     total_undetermined_Update(hObject,handles);
%     go_next_Callback(hObject, eventdata, handles)
% 
% 
% % --- Executes on button press in c_ros.
% function c_ros_Callback(hObject, eventdata, handles)
% 
%     handles.labelName_list{handles.current_idx} = 'bullet';
%     handles.labelId_list(handles.current_idx) = 8;
%     total_flakes_Update(hObject,handles);
%     total_labels_Update(hObject,handles);
%     total_undetermined_Update(hObject,handles);
%     go_next_Callback(hObject, eventdata, handles)
% 
% 
% % --- Executes on button press in c_other.
% function c_other_Callback(hObject, eventdata, handles)
% 
%     handles.labelName_list{handles.current_idx} = 'other';
%     handles.labelId_list(handles.current_idx) = 0;
%     total_flakes_Update(hObject,handles);
%     total_labels_Update(hObject,handles);
%     total_undetermined_Update(hObject,handles);
%     go_next_Callback(hObject, eventdata, handles)

    
 
%% arrow buttons

% --- Executes on button press in go_prev.
function go_prev_Callback(hObject, eventdata, handles)
% hObject    handle to go_prev (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.current_idx > 1
    handles.current_idx = handles.current_idx - 1;
    set(handles.listbox,'Value',handles.current_idx);
    set(handles.is_melting,'Value',handles.isMelting(handles.current_idx));
    load(fullfile(handles.dir_data,handles.file_list{handles.current_idx}));
    handles.current_name = roi.name;
    refresh_image(handles,roi);
    guidata(hObject, handles);
end

% --- Executes on button press in go_next.
function go_next_Callback(hObject, eventdata, handles)
% hObject    handle to go_next (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.current_idx < handles.idx_max
    handles.current_idx = handles.current_idx + 1;
    set(handles.listbox,'Value',handles.current_idx);
    set(handles.is_melting,'Value',handles.isMelting(handles.current_idx));
    load(fullfile(handles.dir_data,handles.file_list{handles.current_idx}));
    handles.current_name = roi.name;
    refresh_image(handles,roi);
    guidata(hObject, handles);
end

% --- Executes on button press in go_end.
function go_end_Callback(hObject, eventdata, handles)
% hObject    handle to go_end (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.current_idx < handles.idx_max
    handles.current_idx = handles.idx_max;
    set(handles.listbox,'Value',handles.current_idx);
    set(handles.is_melting,'Value',handles.isMelting(handles.current_idx));
    load(fullfile(handles.dir_data,handles.file_list{handles.current_idx}));
    handles.current_name = roi.name;
    refresh_image(handles,roi);
    guidata(hObject, handles);    
end

% --- Executes on button press in go_start.
function go_start_Callback(hObject, eventdata, handles)
% hObject    handle to go_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.current_idx > 1
    handles.current_idx = 1;
    set(handles.listbox,'Value',handles.current_idx);
    set(handles.is_melting,'Value',handles.isMelting(handles.current_idx));
    load(fullfile(handles.dir_data,handles.file_list{handles.current_idx}));
    handles.current_name = roi.name;
    refresh_image(handles,roi);
    guidata(hObject,handles);
end

%% listbox functions

%--- Executes during object creation, after setting all properties.
function listbox_CreateFcn(hObject, eventdata, handles)

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in listbox.
function listbox_Callback(hObject, eventdata, handles)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox
get(handles.figure1,'SelectionType');
% If double click
if strcmp(get(handles.figure1,'SelectionType'),'open')
    idx_selected = get(handles.listbox,'Value');
    set(handles.is_melting,'Value',handles.isMelting(handles.current_idx));
    handles.current_idx = idx_selected;
    load(fullfile(handles.dir_data,handles.file_list{handles.current_idx}));
    refresh_image(handles,roi);
    handles.current_name = roi.name;
    guidata(hObject,handles);
 
end

function load_listbox(handles)
set(handles.listbox,'String',handles.file_list);
if ~isempty(handles.file_list)
    set(handles.listbox,'Value',handles.current_idx);
else
    set(handles.listbox,'Value',1);
end


%--- Executes during object creation, after setting all properties.
function text_dir_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
%set(handles.text_dir,'String',handles.dir_data);



% directory button
function disp_current_folder(handles)
set(handles.text_dir,'String',handles.dir_data);



% image button
function refresh_image(handles,roi)
pos =  getpixelposition(handles.axes1);%,'Position');
width = floor(pos(3));
height = floor(pos(4));
if roi.width+2 >= width || roi.height+2 >= height
    imshow(roi.data);
else
    im = uint8(zeros(width,height));
    im(floor(height/2)-floor(roi.height/2)+1:floor(height/2)+ceil(roi.height/2),floor(width/2)-floor(roi.width/2)+1:floor(width/2)+ceil(roi.width/2)) = roi.data;
    imshow(im);
end
if isempty(handles.labelName_list{handles.current_idx})
    title('not labelled yet.');
else
    if handles.isMelting(handles.current_idx)
        title(sprintf('%s / %s / Melting snow',handles.labelName_list{handles.current_idx},handles.rimingName_list{handles.current_idx}),'interpreter','none');
    else    
        title(sprintf('%s / %s',handles.labelName_list{handles.current_idx},handles.rimingName_list{handles.current_idx}),'interpreter','none');
    end
end   



% browse button

% --- Executes on button press in browse_dir.
function browse_dir_Callback(hObject, eventdata, handles)
[~,folder_name] = uigetfile('*.mat','Please click on a file to select the parent folder');
if folder_name ~= 0
    handles.dir_data = folder_name;
    handles.file_list = dir(fullfile(handles.dir_data,'20*.mat'));
    % for Massimo data
    if isempty(handles.file_list)
        handles.file_list = dir(fullfile(handles.dir_data,'ICE*.mat'));
    end
    handles.file_list = {handles.file_list.name}';
    handles.idx_max = length(handles.file_list);
    % randomization
    handles.idx_rnd = randperm(handles.idx_max);
    handles.file_list = handles.file_list(handles.idx_rnd);
    handles.current_idx = 1;
    load(fullfile(handles.dir_data,handles.file_list{handles.current_idx}));
    handles.current_name = roi.name;
    handles.labelId_list = zeros(length(handles.file_list),1);
    handles.labelId_list(handles.labelId_list==0) = -1; 
    handles.labelName_list = cell(length(handles.file_list),1);
    handles.rimingId_list = zeros(length(handles.file_list),1);
    handles.rimingId_list(handles.rimingId_list==0) = -1;
    handles.rimingName_list = cell(length(handles.file_list),1);
    handles.isMelting = false(length(handles.file_list),1);
    set(handles.is_melting,'Value',handles.isMelting(handles.current_idx));
    guidata(hObject,handles);
    refresh_image(handles,roi);
    handles = guidata(hObject);
    load_listbox(handles);
    disp_current_folder(handles);
    total_flakes_Update(hObject,handles);
    total_labels_Update(hObject,handles);
    total_undetermined_Update(hObject,handles);
    console_Update(hObject,handles,sprintf('New folder successfully loaded \n'));
    guidata(hObject,handles);
end
    
%% console functions

% --- Executes during object creation, after setting all properties.
function console_CreateFcn(hObject, eventdata, handles)
% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in console.
function console_Callback(hObject, eventdata, handles)
% Hints: contents = cellstr(get(hObject,'String')) returns console contents as cell array
%        contents{get(hObject,'Value')} returns selected item from console

function console_Init(hObject,handles)
handles.console_list = {'*** MASC Snowflake labelling v.1.1 : updated labelling scheme ***'};
handles.console_line = 1;
set(handles.console,'String',handles.console_list);
set(handles.console,'Value',[]);
guidata(hObject,handles);

function console_Update(hObject,handles,string)
handles.console_list{end+1} = string;
handles.console_line = handles.console_line + 1;
set(handles.console,'String',handles.console_list);
%set(handles.console,'Value',handles.console_line);
guidata(hObject,handles);



%% stats panel functions
function total_flakes_Update(hObject,handles)
set(handles.total_flakes,'String',num2str(length(handles.file_list)));
guidata(hObject,handles);

function total_labels_Update(hObject,handles)
if ~isempty(handles.file_list) > 0
    set(handles.total_labels,'String',sprintf('%u  (%2.1f %%)',sum(handles.labelId_list ~= -1),100*sum(handles.labelId_list ~= -1)/length(handles.file_list)));
else
    set(handles.total_labels,'String',num2str(0));
end
guidata(hObject,handles);

function total_undetermined_Update(hObject,handles)
if ~isempty(handles.file_list) > 0
    set(handles.total_undetermined,'String',sprintf('%u  (%2.1f %%)',sum(handles.labelId_list == 0 & handles.rimingId_list == 0),100*sum(handles.labelId_list == 0 & handles.rimingId_list == 0)/length(handles.file_list)));
else
    set(handles.total_undetermined,'String',num2str(0));
end
guidata(hObject,handles);



function total_flakes_Callback(hObject, eventdata, handles)
% hObject    handle to total_flakes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of total_flakes as text
%        str2double(get(hObject,'String')) returns contents of total_flakes as a double

% --- Executes during object creation, after setting all properties.
function total_flakes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to total_flakes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function total_labels_Callback(hObject, eventdata, handles)
% hObject    handle to total_labels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of total_labels as text
%        str2double(get(hObject,'String')) returns contents of total_labels as a double


% --- Executes during object creation, after setting all properties.
function total_labels_CreateFcn(hObject, eventdata, handles)
% hObject    handle to total_labels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function total_undetermined_Callback(hObject, eventdata, handles)
% hObject    handle to total_undetermined (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of total_undetermined as text
%        str2double(get(hObject,'String')) returns contents of total_undetermined as a double


% --- Executes during object creation, after setting all properties.
function total_undetermined_CreateFcn(hObject, eventdata, handles)
% hObject    handle to total_undetermined (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in save_labels.
function save_labels_Callback(hObject, eventdata, handles)
label.flakes = handles.file_list;
label.className = handles.labelName_list;
label.classId = handles.labelId_list;
label.rimingName = handles.rimingName_list;
label.rimingId = handles.rimingId_list;
label.isMelting = handles.isMelting;
label.idx_rnd = handles.idx_rnd;
if isempty(handles.file_list)
    console_Update(hObject,handles,'Please select a data directory first by clicking on "Browse"');
else
    [fileName,pathName] = uiputfile(fullfile(handles.dir_data,'YourName_labels.mat'),'Save Labels');
    if fileName ~= 0
        save(fullfile(pathName,fileName),'label');
        console_Update(hObject,handles,sprintf('Labels saved successfully in %s',fullfile(pathName,fileName)));
    end
end


% --- Executes on button press in load_labels.
function load_labels_Callback(hObject, eventdata, handles)
if isempty(handles.file_list)
    console_Update(hObject,handles,'Please select a data directory first by clicking on "Browse"');
else    
    [FileName,PathName] = uigetfile('*.mat','Please select a .mat file containing labels');
    if FileName ~= 0
        load(fullfile(PathName,FileName));
        if exist('label') && length(handles.file_list) == length(label.flakes)
            handles.labelName_list = label.className;
            handles.labelId_list = label.classId;
            handles.rimingName_list = label.rimingName;
            handles.rimingId_list = label.rimingId;
            handles.isMelting = label.isMelting;
            handles.idx_rnd = label.idx_rnd;
            handles.file_list = dir(fullfile(handles.dir_data,'20*.mat'));
            % for Massimo data
            if isempty(handles.file_list)
                handles.file_list = dir(fullfile(handles.dir_data,'ICE*.mat'));
            end
            handles.file_list = {handles.file_list.name}';
            handles.file_list = handles.file_list(handles.idx_rnd);
                        
            handles.current_idx = 1;
            set(handles.is_melting,'Value',handles.isMelting(handles.current_idx));
            load(fullfile(handles.dir_data,handles.file_list{handles.current_idx}));
            handles.current_name = roi.name;
            guidata(hObject,handles);
            refresh_image(handles,roi);
            handles = guidata(hObject);
            load_listbox(handles);
            disp_current_folder(handles);

            total_flakes_Update(hObject,handles);
            total_labels_Update(hObject,handles);
            total_undetermined_Update(hObject,handles);
            console_Update(hObject,handles,'Labels loaded successfully. You can now continue to label where you stopped'); 
            guidata(hObject,handles);
        else
            console_Update(hObject,handles,'Please select a valid labels file.');
        end
    end      
end


% --- Executes on key press with focus on figure1 or any of its controls.
function figure1_WindowKeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
switch eventdata.Key
    case 'uparrow'
        go_start_Callback(hObject, eventdata, handles);
    case 'downarrow'
        go_end_Callback(hObject, eventdata, handles);
    case 'leftarrow'
        go_prev_Callback(hObject, eventdata, handles);
    case 'rightarrow'
        go_next_Callback(hObject, eventdata, handles);
    case '1'
        riming_dry_Callback(hObject, eventdata, handles);
    case 'numpad1'
        riming_dry_Callback(hObject, eventdata, handles);
    case '2'
        riming_mrimed_Callback(hObject, eventdata, handles);
    case 'numpad2'
        riming_mrimed_Callback(hObject, eventdata, handles);
    case '3'
        riming_hrimed_Callback(hObject, eventdata, handles);
    case 'numpad3'
        riming_hrimed_Callback(hObject, eventdata, handles);
    case '4'
        riming_graupel_like_Callback(hObject, eventdata, handles);
    case 'numpad4'
        riming_graupel_like_Callback(hObject, eventdata, handles);
    case '5'
        riming_graupel_Callback(hObject, eventdata, handles);
    case 'numpad5'
        riming_graupel_Callback(hObject, eventdata, handles);
    case '0'
        riming_undetermined_Callback(hObject, eventdata, handles);
    case 'numpad0'
        riming_undetermined_Callback(hObject, eventdata, handles);
    case 'add'
        set(handles.is_melting,'Value',1);
        handles.isMelting(handles.current_idx) =1 ;
        guidata(hObject, handles);
    case 'subtract'
        set(handles.is_melting,'Value',0);
        handles.isMelting(handles.current_idx) = 0;
        guidata(hObject, handles);
    otherwise
        disp(eventdata.Key);
end


% --- Executes on button press in class_small.
function class_small_Callback(hObject, eventdata, handles)

    handles.labelName_list{handles.current_idx} = 'small';
    handles.labelId_list(handles.current_idx) = 1;
    if handles.rimingId_list(handles.current_idx) ~= -1
        total_flakes_Update(hObject,handles);
        total_labels_Update(hObject,handles);
        total_undetermined_Update(hObject,handles);
        go_next_Callback(hObject, eventdata, handles);
    else
        guidata(hObject, handles);
    end
    
% --- Executes on button press in class_col.
function class_col_Callback(hObject, eventdata, handles)

    handles.labelName_list{handles.current_idx} = 'column';
    handles.labelId_list(handles.current_idx) = 2;
    if handles.rimingId_list(handles.current_idx) ~= -1
        total_flakes_Update(hObject,handles);
        total_labels_Update(hObject,handles);
        total_undetermined_Update(hObject,handles);
        go_next_Callback(hObject, eventdata, handles);
    else
        guidata(hObject, handles);
    end
    
% --- Executes on button press in class_needle.
function class_needle_Callback(hObject, eventdata, handles)

    handles.labelName_list{handles.current_idx} = 'needle';
    handles.labelId_list(handles.current_idx) = 3;
    if handles.rimingId_list(handles.current_idx) ~= -1
        total_flakes_Update(hObject,handles);
        total_labels_Update(hObject,handles);
        total_undetermined_Update(hObject,handles);
        go_next_Callback(hObject, eventdata, handles);
    else
        guidata(hObject, handles);
    end

% --- Executes on button press in class_plate.
function class_plate_Callback(hObject, eventdata, handles)

    handles.labelName_list{handles.current_idx} = 'plate';
    handles.labelId_list(handles.current_idx) = 4;
    if handles.rimingId_list(handles.current_idx) ~= -1
        total_flakes_Update(hObject,handles);
        total_labels_Update(hObject,handles);
        total_undetermined_Update(hObject,handles);
        go_next_Callback(hObject, eventdata, handles);
    else
        guidata(hObject, handles);
    end

% --- Executes on button press in class_sec_plate.
function class_sec_plate_Callback(hObject, eventdata, handles)

    handles.labelName_list{handles.current_idx} = 'sec_plate';
    handles.labelId_list(handles.current_idx) = 5;
    if handles.rimingId_list(handles.current_idx) ~= -1
        total_flakes_Update(hObject,handles);
        total_labels_Update(hObject,handles);
        total_undetermined_Update(hObject,handles);
        go_next_Callback(hObject, eventdata, handles);
    else
        guidata(hObject, handles);
    end


% --- Executes on button press in class_dendrite.
function class_dendrite_Callback(hObject, eventdata, handles)

    handles.labelName_list{handles.current_idx} = 'dendrite';
    handles.labelId_list(handles.current_idx) = 6;
    if handles.rimingId_list(handles.current_idx) ~= -1
        total_flakes_Update(hObject,handles);
        total_labels_Update(hObject,handles);
        total_undetermined_Update(hObject,handles);
        go_next_Callback(hObject, eventdata, handles);
    else
        guidata(hObject, handles);
    end

% --- Executes on button press in class_agg.
function class_agg_Callback(hObject, eventdata, handles)

    handles.labelName_list{handles.current_idx} = 'aggregate';
    handles.labelId_list(handles.current_idx) = 7;
    if handles.rimingId_list(handles.current_idx) ~= -1
        total_flakes_Update(hObject,handles);
        total_labels_Update(hObject,handles);
        total_undetermined_Update(hObject,handles);
        go_next_Callback(hObject, eventdata, handles);
    else
        guidata(hObject, handles);
    end

% --- Executes on button press in class_combin_col.
function class_combin_col_Callback(hObject, eventdata, handles)

    handles.labelName_list{handles.current_idx} = 'combin_col';
    handles.labelId_list(handles.current_idx) = 8;
    if handles.rimingId_list(handles.current_idx) ~= -1
        total_flakes_Update(hObject,handles);
        total_labels_Update(hObject,handles);
        total_undetermined_Update(hObject,handles);
        go_next_Callback(hObject, eventdata, handles);
    else
        guidata(hObject, handles);
    end

% --- Executes on button press in class_combin_planar.
function class_combin_planar_Callback(hObject, eventdata, handles)

    handles.labelName_list{handles.current_idx} = 'combin_planar';
    handles.labelId_list(handles.current_idx) = 9;
    if handles.rimingId_list(handles.current_idx) ~= -1
        total_flakes_Update(hObject,handles);
        total_labels_Update(hObject,handles);
        total_undetermined_Update(hObject,handles);
        go_next_Callback(hObject, eventdata, handles);
    else
        guidata(hObject, handles);
    end

% --- Executes on button press in class_combin_mixte.
function class_combin_mixte_Callback(hObject, eventdata, handles)

    handles.labelName_list{handles.current_idx} = 'combin_mixte';
    handles.labelId_list(handles.current_idx) = 10;
    if handles.rimingId_list(handles.current_idx) ~= -1
        total_flakes_Update(hObject,handles);
        total_labels_Update(hObject,handles);
        total_undetermined_Update(hObject,handles);
        go_next_Callback(hObject, eventdata, handles);
    else
        guidata(hObject, handles);
    end
    
% --- Executes on button press in class_graupel.
function class_graupel_Callback(hObject, eventdata, handles)

    handles.labelName_list{handles.current_idx} = 'graupel';
    handles.labelId_list(handles.current_idx) = 11;
    handles.rimingName_list{handles.current_idx} = 'graupel';
    handles.rimingId_list(handles.current_idx) = 5;
    if handles.rimingId_list(handles.current_idx) ~= -1
        total_flakes_Update(hObject,handles);
        total_labels_Update(hObject,handles);
        total_undetermined_Update(hObject,handles);
        go_next_Callback(hObject, eventdata, handles);
    else
        guidata(hObject, handles);
    end

% --- Executes on button press in class_undetermined.
function class_undetermined_Callback(hObject, eventdata, handles)

    handles.labelName_list{handles.current_idx} = 'undetermined';
    handles.labelId_list(handles.current_idx) = 0;
    if handles.rimingId_list(handles.current_idx) ~= -1
        total_flakes_Update(hObject,handles);
        total_labels_Update(hObject,handles);
        total_undetermined_Update(hObject,handles);
        go_next_Callback(hObject, eventdata, handles);
    else
        guidata(hObject, handles);
    end

% --- Executes on button press in riming_dry.
function riming_dry_Callback(hObject, eventdata, handles)

    handles.rimingName_list{handles.current_idx} = 'none';
    handles.rimingId_list(handles.current_idx) = 1;
    if handles.labelId_list(handles.current_idx) ~= -1
        total_flakes_Update(hObject,handles);
        total_labels_Update(hObject,handles);
        total_undetermined_Update(hObject,handles);
        go_next_Callback(hObject, eventdata, handles);
    else
        guidata(hObject, handles);
    end

% --- Executes on button press in riming_mrimed.
function riming_mrimed_Callback(hObject, eventdata, handles)

    handles.rimingName_list{handles.current_idx} = 'm_rimed';
    handles.rimingId_list(handles.current_idx) = 2;
    if handles.labelId_list(handles.current_idx) ~= -1
        total_flakes_Update(hObject,handles);
        total_labels_Update(hObject,handles);
        total_undetermined_Update(hObject,handles);
        go_next_Callback(hObject, eventdata, handles);
    else
        guidata(hObject, handles);
    end

% --- Executes on button press in riming_hrimed.
function riming_hrimed_Callback(hObject, eventdata, handles)

    handles.rimingName_list{handles.current_idx} = 'h_rimed';
    handles.rimingId_list(handles.current_idx) = 3;
    if handles.labelId_list(handles.current_idx) ~= -1
        total_flakes_Update(hObject,handles);
        total_labels_Update(hObject,handles);
        total_undetermined_Update(hObject,handles);
        go_next_Callback(hObject, eventdata, handles);
    else
        guidata(hObject, handles);
    end

% --- Executes on button press in riming_graupel_like.
function riming_graupel_like_Callback(hObject, eventdata, handles)

    handles.rimingName_list{handles.current_idx} = 'graupel_like';
    handles.rimingId_list(handles.current_idx) = 4;
    if handles.labelId_list(handles.current_idx) ~= -1
        total_flakes_Update(hObject,handles);
        total_labels_Update(hObject,handles);
        total_undetermined_Update(hObject,handles);
        go_next_Callback(hObject, eventdata, handles);
    else
        guidata(hObject, handles);
    end

% --- Executes on button press in riming_graupel.
function riming_graupel_Callback(hObject, eventdata, handles)

    handles.rimingName_list{handles.current_idx} = 'graupel';
    handles.rimingId_list(handles.current_idx) = 5;
    handles.labelName_list{handles.current_idx} = 'graupel';
    handles.labelId_list(handles.current_idx) = 11;
    if handles.labelId_list(handles.current_idx) ~= -1
        total_flakes_Update(hObject,handles);
        total_labels_Update(hObject,handles);
        total_undetermined_Update(hObject,handles);
        go_next_Callback(hObject, eventdata, handles);
    else
        guidata(hObject, handles);
    end

% --- Executes on button press in riming_undetermined.
function riming_undetermined_Callback(hObject, eventdata, handles)

    handles.rimingName_list{handles.current_idx} = 'undetermined';
    handles.rimingId_list(handles.current_idx) = 0;
    if handles.labelId_list(handles.current_idx) ~= -1
        total_flakes_Update(hObject,handles);
        total_labels_Update(hObject,handles);
        total_undetermined_Update(hObject,handles);
        go_next_Callback(hObject, eventdata, handles);
    else
        guidata(hObject, handles);
    end




% --- Executes on button press in is_melting.
function is_melting_Callback(hObject, eventdata, handles)

    handles.isMelting(handles.current_idx) = get(hObject,'Value');
    guidata(hObject, handles);
    

% Hint: get(hObject,'Value') returns toggle state of is_melting

function varargout = megui(varargin)
% MEGUI MATLAB code for megui.fig
%      MEGUI, by itself, creates a new MEGUI or raises the existing
%      singleton*.
%
%      H = MEGUI returns the handle to a new MEGUI or the handle to
%      the existing singleton*.
%
%      MEGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MEGUI.M with the given input arguments.
%
%      MEGUI('Property','Value',...) creates a new MEGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before megui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to megui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help megui

% Last Modified by GUIDE v2.5 24-Dec-2024 19:50:26

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @megui_OpeningFcn, ...
                   'gui_OutputFcn',  @megui_OutputFcn, ...
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


% --- Executes just before megui is made visible.
function megui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to megui (see VARARGIN)

% Choose default command line output for megui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes megui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = megui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --------------------------------------------------------------------
function File_Callback(hObject, eventdata, handles)
% hObject    handle to File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function transformation_Callback(hObject, eventdata, handles)
% hObject    handle to transformation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function enhancement_Callback(hObject, eventdata, handles)
% hObject    handle to enhancement (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function tuxiangfenge_Callback(hObject, eventdata, handles)
% hObject    handle to tuxiangfenge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function tuxianglvbo_Callback(hObject, eventdata, handles)
% hObject    handle to tuxianglvbo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function edge_Callback(hObject, eventdata, handles)
% hObject    handle to edge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Sobel_Callback(hObject, eventdata, handles)
% hObject    handle to Sobel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global im;
global out;

rgb = im;
mysize = size(rgb);
if numel(mysize)>2
  im1 = rgb2gray(rgb); %将彩色图像转换为灰度图像
else
im1 = rgb;
end

out = edge(im1,'Sobel',0.06);
axes(handles.axes2); 
imshow(out);title('Sobel算子边缘检测');

% --------------------------------------------------------------------
function Canny_Callback(hObject, eventdata, handles)
% hObject    handle to Canny (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global im;
global out;

rgb = im;
mysize = size(rgb);
if numel(mysize)>2
  im1 = rgb2gray(rgb); %将彩色图像转换为灰度图像
else
im1 = rgb;
end

out = edge(im1,'Canny',0.06);
axes(handles.axes2); 
imshow(out);title('Canny算子边缘检测');

% --------------------------------------------------------------------
function Log_Callback(hObject, eventdata, handles)
% hObject    handle to Log (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global im;
global out;

rgb = im;
mysize = size(rgb);
if numel(mysize)>2
  im1 = rgb2gray(rgb); %将彩色图像转换为灰度图像
else
im1 = rgb;
end

out = edge(im1,'Log',0.012);
axes(handles.axes2); 
imshow(out);title('Log算子边缘检测');

% --------------------------------------------------------------------
function Prewitt_Callback(hObject, eventdata, handles)
% hObject    handle to Prewitt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global im;
global out;

rgb = im;
mysize = size(rgb);
if numel(mysize)>2
  im1 = rgb2gray(rgb); %将彩色图像转换为灰度图像
else
im1 = rgb;
end

out = edge(im1,'Prewitt',0.06);
axes(handles.axes2); 
imshow(out);title('Prewitt算子边缘检测');

% --------------------------------------------------------------------
function Roberts_Callback(hObject, eventdata, handles)
% hObject    handle to Roberts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global im;
global out;

rgb = im;
mysize = size(rgb);
if numel(mysize)>2
  im1 = rgb2gray(rgb); %将彩色图像转换为灰度图像
else
im1 = rgb;
end

out = edge(im1,'Roberts',0.06);
axes(handles.axes2); 
imshow(out);title('Roberts算子边缘检测');


% --------------------------------------------------------------------
function zhongzhilvbo_Callback(hObject, eventdata, handles)
% hObject    handle to zhongzhilvbo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global im;
global out;

rgb = im;
mysize = size(rgb);
if numel(mysize)>2
  im = rgb2gray(rgb); %将彩色图像转换为灰度图像
else
im = rgb;
end
fn=imnoise(im,'salt & pepper',0.2);%用函数imnosie产生椒盐噪声，0.2代表图中白点黑点出现的概率为0.2
out = medfilt2(fn);%中值滤波
axes(handles.axes2); 
imshow(out);title('中值滤波处理后的图像'); 

% --------------------------------------------------------------------
function batewosi_Callback(hObject, eventdata, handles)
% hObject    handle to batewosi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global im;
global out;

rgb = im;
mysize = size(rgb);
if numel(mysize)>2
  I1 = rgb2gray(rgb); %将彩色图像转换为灰度图像
else
I1 = rgb;
end
m = double(I1);
f = fft2(m);
f = fftshift(f);
[N1,N2] = size(f);
n1 = round(N1/2);
n2 = round(N2/2);
n = 2;
d0 = 10;
for i = 1:N1
    for j = 1:N2
        d = sqrt((i-n1)^2+(j-n2)^2);
        h = (1/(1+(d0/d)^(2*n)))+0.5;
        y(i,j) = h*f(i,j);
    end
end
y = ifftshift(y);
A = ifft2(y);
out = uint8(real(A));
axes(handles.axes2); 
imshow(out);title('巴特沃斯高通滤波后的图像'); 

% --------------------------------------------------------------------
function dian_Callback(hObject, eventdata, handles)
% hObject    handle to dian (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global im;
global out;

rgb = im;
mysize = size(rgb);
if numel(mysize)>2
  im1 = rgb2gray(rgb); %将彩色图像转换为灰度图像
else
im1 = rgb;
end

w = [-1 -1 -1;-1 8 -1;-1 -1 -1];%给定模板
out = ordfilt2(im1,5*5,ones(5,5))-ordfilt2(im1,1,ones(5,5));%在这里采用了5*5模板进行差值处理
T2 = max(out(:));%同理
out = out>=T2/2;%同理
axes(handles.axes2); 
imshow(out);title('点检测');

% --------------------------------------------------------------------
function xian_Callback(hObject, eventdata, handles)
% hObject    handle to xian (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global im;
global out;

rgb = im;
mysize = size(rgb);
if numel(mysize)>2
  im1 = rgb2gray(rgb); %将彩色图像转换为灰度图像
else
im1 = rgb;
end

w = [2 -1 -1;-1 2 -1;-1 -1 2];          % -45°方向检测线
out = imfilter(double(im1),w);
axes(handles.axes2); 
imshow(out);title('线检测');

% --------------------------------------------------------------------
function speckle_Callback(hObject, eventdata, handles)
% hObject    handle to speckle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function huiduhua_Callback(hObject, eventdata, handles)
% hObject    handle to huiduhua (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global im
global BW  % 定义全局变量
axes(handles.axes2);

% 检查 im 是否为 RGB 图像
if size(im, 3) == 3
    BW = rgb2gray(im);  % 如果是 RGB 图像，则转换为灰度图像
else
    BW = im;  % 如果已经是灰度图像，则直接使用
end

imshow(BW);  % 显示灰度图像
% --------------------------------------------------------------------
function pinghualvbo_Callback(hObject, eventdata, handles)
% hObject    handle to pinghualvbo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global im;
global out;

Inoised = imnoise(im,'gaussian',0.1,0.005);%对图像进行高斯噪声加噪
%制定卷积核
h=ones(3,3)/5;
h(1,1) = 0;
h(1,3) = 0;
h(3,1) = 0;
h(1,3) = 0;
%平滑运算
out = imfilter(Inoised,h);
axes(handles.axes2); 
imshow(out);title('平滑滤波处理后的图像'); 

% --------------------------------------------------------------------
function ruihua_Callback(hObject, eventdata, handles)
% hObject    handle to ruihua (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global im;
global out;

rgb = im;
mysize = size(rgb);
if numel(mysize)>2
  I1 = rgb2gray(rgb); %将彩色图像转换为灰度图像
else
I1 = rgb;
end
model=[-1,0,1;
       -2,0,2;
       -1,0,1];
[m,n]=size(I1);
I2=double(I1);

for i=2:m-1
    for j=2:n-1
        I2(i,j)=I1(i+1,j+1)+2*I1(i+1,j)+I1(i+1,j-1)-I1(i-1,j+1)-2*I1(i-1,j)-I1(i-1,j-1);
    end
end
I2 = I2 + double(I1);
out = I2;
axes(handles.axes2); 
imshow(uint8(out));title('锐化后的图像'); 

% --------------------------------------------------------------------
function xuanzhuan_Callback(hObject, eventdata, handles)
% hObject    handle to xuanzhuan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global im
global out  %定义全局变量 

% A=getimage(handles.axes1);
A  = im;
axes(handles.axes2); 
prompt = {'度数：'};
def={'90'};
answer = inputdlg(prompt,'请输入：',1,def);
if ~isempty(answer)
a = str2num(answer{1});
J = imrotate(A,360-a);
out = J;
axes(handles.axes2); 
imshow(out,[]);title('旋转后图像'); 
end

% --------------------------------------------------------------------
function DFT_Callback(hObject, eventdata, handles)
% hObject    handle to DFT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global im
global out  

rgb = im;
mysize = size(rgb);
if numel(mysize)>2
  im = rgb2gray(rgb); %将彩色图像转换为灰度图像
else
im = rgb;
end

I1 = im2double(im);
I2 = fft2(I1);
I3 = fftshift(I2);
out = log(abs(I3)+1); 
axes(handles.axes2); 
imshow(out,[]);title('离散傅里叶变换'); 

% --------------------------------------------------------------------
function Open_Callback(hObject, eventdata, handles)
% hObject    handle to Open (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global im   %定义一个全局变量im
global im2
[filename,pathname]=...
    uigetfile({'*.*';'*.bmp';'*.tif';'*.png'},'select picture');  %选择图片路径
str=[pathname filename];  %合成路径+文件名
im=imread(str);   %读取图片
im2=im;
axes(handles.axes1);  %使用第一个axes
imshow(im);  %显示图片

% --------------------------------------------------------------------
function Save_Callback(hObject, eventdata, handles)
% hObject    handle to Save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global BW 
set(handles.axes2,'HandleVisibility','ON');
axes(handles.axes2);
[filename,pathname]=uiputfile({'*.jpg';'*.bmp';'*.tif';'*.*'},'save image as');
file=strcat(pathname,filename);
BW=getimage(gca);
imwrite(BW,file);
set(handles.axes2,'HandleVisibility','Off');

% --------------------------------------------------------------------
function Quit_Callback(hObject, eventdata, handles)
% hObject    handle to Quit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(gcf)  %关闭当前Figure窗口句柄


% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
set(gca,'XColor',get(gca,'Color')) ;% 这两行代码功能：将坐标轴和坐标刻度转为白色
set(gca,'YColor',get(gca,'Color'));
 
set(gca,'XTickLabel',[]); % 这两行代码功能：去除坐标刻度
set(gca,'YTickLabel',[]);
% Hint: place code in OpeningFcn to populate axes1


% --- Executes during object creation, after setting all properties.
function axes2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
set(gca,'XColor',get(gca,'Color')) ;% 这两行代码功能：将坐标轴和坐标刻度转为白色
set(gca,'YColor',get(gca,'Color'));
 
set(gca,'XTickLabel',[]); % 这两行代码功能：去除坐标刻度
set(gca,'YTickLabel',[]);
% Hint: place code in OpeningFcn to populate axes2


% --------------------------------------------------------------------
function process_Callback(hObject, eventdata, handles)
% hObject    handle to process (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function huiduzhifangtu_Callback(hObject, eventdata, handles)
% hObject    handle to huiduzhifangtu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 global im;
    grayImage = im2gray(im); % 将图像转换为灰度（兼容 RGB 和灰度图像）
    axes(handles.axes2); % 显示在第二个 axes
    imhist(grayImage); % 绘制灰度直方图
    title('灰度直方图');
% --------------------------------------------------------------------
function zftjunhenghua_Callback(hObject, eventdata, handles)
% hObject    handle to zftjunhenghua (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    global im;
    grayImage = im2gray(im); % 使用 im2gray 以确保图像为灰度
    axes(handles.axes2); % 显示在第二个 axes
    eqImage = histeq(grayImage); % 进行直方图均衡化
    imshow(eqImage);
    title('直方图均衡化后的图像');

% --------------------------------------------------------------------
function zftpipei_Callback(hObject, eventdata, handles)
% hObject    handle to zftpipei (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    global im;
    [filename, pathname] = uigetfile({'*.jpg;*.png;*.bmp'}, '选择参考图像');
    if isequal(filename, 0)
        return;
    end
    refImage = imread(fullfile(pathname, filename)); % 读取参考图像
    grayImage = im2gray(im); % 将待处理图像转换为灰度
    refGrayImage = im2gray(refImage); % 将参考图像转换为灰度
    matchedImage = imhistmatch(grayImage, refGrayImage); % 直方图匹配
    axes(handles.axes2); % 显示在第二个 axes
    imshow(matchedImage);
    title('直方图匹配后的图像');


% --------------------------------------------------------------------
function noise_Callback(hObject, eventdata, handles)
% hObject    handle to noise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function gosi_Callback(hObject, eventdata, handles)
% hObject    handle to gosi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global im
global out  %定义全局变量 

I = im2double(im);
% 弹出输入框，让用户输入噪声大小
prompt = {'请输入噪声大小 (0 到 1)：'};
def = {'0.02'};  % 默认噪声大小
answer = inputdlg(prompt, '输入噪声大小', 1, def);

if ~isempty(answer)
    noiseLevel = str2double(answer{1});  % 获取用户输入的噪声大小
    if isnan(noiseLevel) || noiseLevel < 0 || noiseLevel > 1
        msgbox('噪声大小必须是 0 到 1 之间的数值。', '错误', 'error');
        return;
    end

J = imnoise(I,'gaussian',noiseLevel);
out = J;
axes(handles.axes2); 
imshow(out);title('添加高斯噪声后的图像'); 
end
% --------------------------------------------------------------------
function bosong_Callback(hObject, eventdata, handles)
% hObject    handle to bosong (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global im
global out  %定义全局变量 

I = im2double(im);
% 弹出输入框，让用户输入噪声大小
prompt = {'请输入噪声大小 (0 到 1)：'};
def = {'0.02'};  % 默认噪声大小
answer = inputdlg(prompt, '输入噪声大小', 1, def);

if ~isempty(answer)
    noiseLevel = str2double(answer{1});  % 获取用户输入的噪声大小
    if isnan(noiseLevel) || noiseLevel < 0 || noiseLevel > 1
        msgbox('噪声大小必须是 0 到 1 之间的数值。', '错误', 'error');
        return;
    end
 %调整图像强度并添加泊松噪声
I_scaled = I * noiseLevel;        % 缩放图像
J = imnoise(I_scaled,'poisson');
% J = J  / noiseLevel;
out = J;
axes(handles.axes2); 
imshow(out);title('添加泊松噪声后的图像'); 
end
% --------------------------------------------------------------------
function yan_Callback(hObject, eventdata, handles)
% hObject    handle to yan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global im
global out  %定义全局变量 

I = im2double(im);
% 弹出输入框，让用户输入噪声大小
prompt = {'请输入噪声大小 (0 到 1)：'};
def = {'0.02'};  % 默认噪声大小
answer = inputdlg(prompt, '输入噪声大小', 1, def);

if ~isempty(answer)
    noiseLevel = str2double(answer{1});  % 获取用户输入的噪声大小
    if isnan(noiseLevel) || noiseLevel < 0 || noiseLevel > 1
        msgbox('噪声大小必须是 0 到 1 之间的数值。', '错误', 'error');
        return;
    end

J = imnoise(I,'salt',noiseLevel);
out = J;
axes(handles.axes2); 
imshow(out);title('添加盐噪声后的图像'); 
end
% --------------------------------------------------------------------
function jiao_Callback(hObject, eventdata, handles)
% hObject    handle to jiao (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global im
global out  %定义全局变量 

I = im2double(im);
% 弹出输入框，让用户输入噪声大小
prompt = {'请输入噪声大小 (0 到 1)：'};
def = {'0.02'};  % 默认噪声大小
answer = inputdlg(prompt, '输入噪声大小', 1, def);

if ~isempty(answer)
    noiseLevel = str2double(answer{1});  % 获取用户输入的噪声大小
    if isnan(noiseLevel) || noiseLevel < 0 || noiseLevel > 1
        msgbox('噪声大小必须是 0 到 1 之间的数值。', '错误', 'error');
        return;
    end

J = imnoise(I,'speckle',noiseLevel);
out = J;
axes(handles.axes2); 
imshow(out);title('添加胡椒噪声后的图像');
end


% --------------------------------------------------------------------
function tztq_Callback(hObject, eventdata, handles)
% hObject    handle to tztq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function tqmb_Callback(hObject, eventdata, handles)
% hObject    handle to tqmb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    global im;  % 假设图像已经在 global 变量 im 中加载
    
    % 确保图像已加载
    if isempty(im)
        msgbox('请先加载图像！', '错误', 'error');
        return;
    end
    
    % 获取原始图像
    img = im;
    
    % 将彩色图像转换为灰度图像（如果是彩色图像）
    if size(img, 3) == 3
        img_gray = rgb2gray(img);  % 将彩色图像转换为灰度图像
    else
        img_gray = img;  % 如果已经是灰度图像，直接使用
    end
    
    % 使用 Otsu 阈值方法自动计算最佳阈值
    threshold = graythresh(img_gray);
    
    % 对图像进行二值化处理
    img_binary = imbinarize(img_gray, threshold);
    
    % 提取目标区域
    target_image = img;  % 默认情况下，目标图像就是原始图像
    if size(img, 3) == 3
        % 如果是彩色图像，将背景置为黑色（掩膜外的部分设为黑色）
        target_image(repmat(~img_binary, [1, 1, 3])) = 0;
    else
        % 如果是灰度图像，将背景置为黑色
        target_image(~img_binary) = 0;
    end
    
    % 将提取的目标图像保存到 handles 结构体中，以便后续使用
    handles.target_image = target_image;
    guidata(hObject, handles);  % 更新 handles 结构体
    
    % 在 GUI 中显示提取的目标图像
    axes(handles.axes2);  % 设置显示区域
    imshow(target_image);  % 显示目标图像
    title('提取的目标图像');  % 设置标题


% --------------------------------------------------------------------
function lbp_Callback(hObject, eventdata, handles)
% hObject    handle to lbp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    global im;
    
    % 检查图像是否已经加载
    if isempty(im)
        msgbox('请先加载图像！', '错误', 'error');
        return;
    end
    
    % 提取原始图像的 LBP 特征并显示到文本框
    extract_LBP_features();

function extract_LBP_features()
    global im;
    
    % 如果图像是彩色图像，转换为灰度图像
    if size(im, 3) == 3
        img_gray = rgb2gray(im);  % 转换为灰度图
    else
        img_gray = im;  % 如果已经是灰度图像，直接使用
    end
    
    % 提取 LBP 特征
    lbp_features = extractLBPFeatures(img_gray, 'CellSize', [32, 32]);
    
    % 将 LBP 特征转换为字符串，以便显示在文本框中
    lbp_str = num2str(lbp_features);  % 将特征向量转换为字符串
    
    disp(lbp_str);



% --------------------------------------------------------------------
function hog_Callback(hObject, eventdata, handles)
% hObject    handle to hog (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    global im;
    
    % 检查图像是否已经加载
    if isempty(im)
        msgbox('请先加载图像！', '错误', 'error');
        return;
    end
    
    % 提取原始图像的 HOG 特征并显示到文本框
    extract_HOG_features();

function extract_HOG_features()
    global im;
    
    % 如果图像是彩色图像，转换为灰度图像
    if size(im, 3) == 3
        img_gray = rgb2gray(im);  % 转换为灰度图
    else
        img_gray = im;  % 如果已经是灰度图像，直接使用
    end
    
    % 提取 HOG 特征
    [hog_features, visualization] = extractHOGFeatures(img_gray, 'CellSize', [8, 8]);
    
    % 将 HOG 特征转换为字符串，以便显示在文本框中
    hog_str = num2str(hog_features);  % 将特征向量转换为字符串
    
    % 显示 HOG 特征到文本框
    disp(hog_str);

% --------------------------------------------------------------------
function xxbh_Callback(hObject, eventdata, handles)
% hObject    handle to xxbh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    global im;  % 使用全局变量 im
    % 如果是彩色图像，转换为灰度图像
    if size(im, 3) == 3
        im_gray = rgb2gray(im);  % 将彩色图像转换为灰度图像
    else
        im_gray = im;  % 如果已经是灰度图，直接使用
    end
  
    % 获取图像的最小值和最大值
    min_val = double(min(im_gray(:)));
    max_val = double(max(im_gray(:)));
    
    % 线性变换：将灰度值映射到 [0, 255] 范围
    a = 0;  % 目标范围的下界
    b = 255;  % 目标范围的上界
    linear_enhanced = uint8((double(im_gray) - min_val) / (max_val - min_val) * (b - a) + a);
    
    % 显示线性变换后的图像
    axes(handles.axes2);  % 设置显示区域
    imshow(linear_enhanced);  % 显示线性变换后的图像
    title('线性对比度增强');
    % 将处理后的图像保存到 handles 中，以便后续使用
    handles.im_gray = im_gray;
    handles.linear_enhanced = linear_enhanced;
    
    % 更新 handles 结构体
    guidata(hObject, handles);

% --------------------------------------------------------------------
function zsbh_Callback(hObject, eventdata, handles)
% hObject    handle to zsbh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    global im;  % 使用全局变量 im
    % 如果是彩色图像，转换为灰度图像
    if size(im, 3) == 3
        im_gray = rgb2gray(im);  % 将彩色图像转换为灰度图像
    else
        im_gray = im;  % 如果已经是灰度图，直接使用
    end
    
    % 将图像数据转换为 double 类型，方便进行运算
    im_gray = double(im_gray);
    
    % 设置指数变换的参数（可以根据需要调整）
    gamma = 0.5;  % 设置指数变换的 γ 值，通常 0 < γ < 1 或 γ > 1
    c = 1;  % 常数，通常设为 1
    
    % 执行指数变换
    im_exp = c * (im_gray .^ gamma);
    
    % 归一化到 [0, 255] 范围，并转换为 uint8 类型
    im_exp = uint8(255 * mat2gray(im_exp));
    
    % 显示指数变换后的图像
    axes(handles.axes2);  % 设置显示区域
    imshow(im_exp);  % 显示指数变换后的图像
    title('指数变换后的图像');
    
    % 将处理后的图像保存到 handles 中，以便后续使用
    handles.im_gray = im_gray;
    handles.im_exp = im_exp;
    
    % 更新 handles 结构体
    guidata(hObject, handles);

% --------------------------------------------------------------------
function dsbh_Callback(hObject, eventdata, handles)
% hObject    handle to dsbh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    global im;  % 使用全局变量 im
    % 如果是彩色图像，转换为灰度图像
    if size(im, 3) == 3
        im_gray = rgb2gray(im);  % 将彩色图像转换为灰度图像
    else
        im_gray = im;  % 如果已经是灰度图，直接使用
    end
    
    % 将图像数据转换为 double 类型，方便进行运算
    im_gray = double(im_gray);
    
    % 设置对数变换的参数（可以根据需要调整）
    c = 1;  % 常数，通常设为 1，控制图像对比度
    % 执行对数变换
    im_log = c * log(1 + im_gray);
    
    % 归一化到 [0, 255] 范围，并转换为 uint8 类型
    im_log = uint8(255 * mat2gray(im_log));
    
    % 显示对数变换后的图像
    axes(handles.axes2);  % 设置显示区域
    imshow(im_log);  % 显示对数变换后的图像
    title('对数变换后的图像');
    
    % 将处理后的图像保存到 handles 中，以便后续使用
    handles.im_gray = im_gray;
    handles.im_log = im_log;
    
    % 更新 handles 结构体
    guidata(hObject, handles);

% --------------------------------------------------------------------
function sx_Callback(hObject, eventdata, handles)
% hObject    handle to sx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% 检查是否已加载图像
    global im im2;
    % 如果 im 为空，提示用户加载图像
    if isempty(im)
        msgbox('请先加载图片！', '错误', 'error');
        return;
    end
    % 设置缩小的比例（如缩小为原图的50%）
    scaleFactor = 0.1;
    % 使用 imresize 函数缩小图像
    im_resized = imresize(im, scaleFactor);
    % 将缩小后的图像保存到 im2（如果需要）
    im2 = im_resized;
    % 显示缩小后的图像
    axes(handles.axes2); % 设置显示区域为 axes1
    imshow(im_resized);  % 显示缩小后的图像
    title('缩小后的图像'); % 添加标题

% --------------------------------------------------------------------
function fd_Callback(hObject, eventdata, handles)
% hObject    handle to fd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    global im;
    
    % 检查是否已加载图像
    if isempty(im)
        msgbox('请先加载图片！', '错误', 'error');
        return;
    end
    
    % 放大图像显示，设置适当的缩放比例
    magnification = 800;  % 设置显示图像的放大比例，200 表示 200%
    
    % 在指定的 axes 上显示图像并放大显示
    axes(handles.axes2);  % 设置显示区域为 axes1
    imshow(im, 'InitialMagnification', magnification);   % 设置放大显示
    title('放大后的图像显示');



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


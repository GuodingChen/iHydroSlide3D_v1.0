classdef Variable_initialize
    properties
        % Hydrological simulation time control
        startYear;startMonth;startDay;startHour;startMinute;startSecond;
        endYear;endMonth;endDay;endHour;endMinute;endSecond;
        % Choose the model core
        modelCore;
        % Forcing data: rain and pet
        rainPath;petPath;rainFormat;petFormat;
        % basic Hydrological input
        Hydro_basicPath; Hydro_basicFormat;
        % related parameters
        paramPath;paramFormat; icsPath;icsFormat; obsPath;
        % output and state file 
        outPath; statePath;stateFormat; resultPath;resultFormat;
        % soil downscaling data
        DownscalingBasicPath; DownscalingBasicFormat;
        % landslide input
        Landslide_basicFormat; Landslide_basicPath;
        UseMaxLocalNcores; UserDefinedNcores;
        Landslide_density; Divide_tile_number;
        min_ae; max_ae;
        min_be; max_be;
        OutputEachStep;
        LandslideOutStart_data; LandslideOutEnd_data;
        timeStep; %outFreq;
        saveState;loadState;
        
        
        RainFactType;KsatType; WMType; BType; IMType; KEType; coeMType;
        %expMType; coeRType; coeSType; KSType; KIType;
        RainFact;Ksat; WM; B; IM; KE; coeM;
        %expM; coeR; coeS; KS; KI;
        
        W0Type;  RType; RSType; %SS0Type; 
        SI0Type; W0;  R; RS;SI0;  %SS0; 
        
        NOutPixs; OutPixNames; OutPixLats; OutPixLons; OutPixAreas;
        
        output_runoff; output_W;output_SM; %output_ExcS;output_ExcI;output_pet;
        output_FS;output_rain;output_epot;output_eact;output_rs;
        output_ri;output_infil;output_tot_infil;
        output_FS3D;
        output_PF;
        output_Volume;
        output_Area;
        
    end
    
    methods
        function obj = Variable_initialize(gFile)
           gfileID = fopen(gFile);
           commentSymbol='#';
           % readLine is a new function defined in this classdef
           % -------------------read the hydrological input--------------
           % read the start and end time of the hydrological simulation
           % the file should be read in order with the controlfile
           obj.startYear = Variable_initialize.readLine(gfileID,'StartYear',commentSymbol,'double');
           obj.startMonth=Variable_initialize.readLine(gfileID,'StartMonth',commentSymbol,'double');
           obj.startDay=Variable_initialize.readLine(gfileID,'StartDay',commentSymbol,'double');
           obj.startHour=Variable_initialize.readLine(gfileID,'StartHour',commentSymbol,'double');
           obj.startMinute=Variable_initialize.readLine(gfileID,'StartMinute',commentSymbol,'double');
           obj.startSecond=Variable_initialize.readLine(gfileID,'StartSecond',commentSymbol,'double');
           
           obj.endYear=Variable_initialize.readLine(gfileID,'EndYear',commentSymbol,'double');
           obj.endMonth=Variable_initialize.readLine(gfileID,'EndMonth',commentSymbol,'double');
           obj.endDay=Variable_initialize.readLine(gfileID,'EndDay',commentSymbol,'double');
           obj.endHour=Variable_initialize.readLine(gfileID,'EndHour',commentSymbol,'double');
           obj.endMinute=Variable_initialize.readLine(gfileID,'EndMinute',commentSymbol,'double');
           obj.endSecond=Variable_initialize.readLine(gfileID,'EndSecond',commentSymbol,'double');
           
           % define the timestep interval
           obj.timeStep=Variable_initialize.readLine(gfileID,'TimeStep',commentSymbol,'double');
%            obj.outFreq=Variable_initialize.readLine(gfileID,'OutFreq',commentSymbol,'double');
  
           obj.modelCore=Variable_initialize.readLine(gfileID,'ModelCore',commentSymbol,'string');
           obj.saveState=Variable_initialize.readLine(gfileID,'SaveState',commentSymbol,'boolean');
           obj.loadState=Variable_initialize.readLine(gfileID,'LoadState',commentSymbol,'boolean');
           
           obj.W0Type=Variable_initialize.readLine(gfileID,'W0Type',commentSymbol,'string');
           obj.W0=Variable_initialize.readLine(gfileID,'W0',commentSymbol,'double');
           obj.RType=Variable_initialize.readLine(gfileID,'RType',commentSymbol,'string');
           obj.R=Variable_initialize.readLine(gfileID,'R',commentSymbol,'double');
           obj.RSType=Variable_initialize.readLine(gfileID,'RSType',commentSymbol,'string');
           obj.RS=Variable_initialize.readLine(gfileID,'RS',commentSymbol,'double');
           %obj.SS0Type=Variable_initialize.readLine(gfileID,'SS0Type',commentSymbol,'string');
           %obj.SS0=Variable_initialize.readLine(gfileID,'SS0',commentSymbol,'double');
           obj.SI0Type=Variable_initialize.readLine(gfileID,'SI0Type',commentSymbol,'string');
           obj.SI0=Variable_initialize.readLine(gfileID,'SI0',commentSymbol,'double');
           
           
            
           obj.Hydro_basicFormat=Variable_initialize.readLine(gfileID,'HydroBasicFormat',commentSymbol,'string');
           obj.Hydro_basicPath=Variable_initialize.readLine(gfileID,'HydroBasicPath',commentSymbol,'string');
           obj.paramFormat=Variable_initialize.readLine(gfileID,'ParamFormat',commentSymbol,'string');
           obj.paramPath=Variable_initialize.readLine(gfileID,'ParamPath',commentSymbol,'string');
           obj.stateFormat=Variable_initialize.readLine(gfileID,'StateFormat',commentSymbol,'string');
           obj.statePath=Variable_initialize.readLine(gfileID,'StatePath',commentSymbol,'string');
           obj.icsFormat=Variable_initialize.readLine(gfileID,'ICSFormat',commentSymbol,'string');
           obj.icsPath=Variable_initialize.readLine(gfileID,'ICSPath',commentSymbol,'string');
           obj.rainFormat=Variable_initialize.readLine(gfileID,'RainFormat',commentSymbol,'string');
           obj.rainPath=Variable_initialize.readLine(gfileID,'RainPath',commentSymbol,'string');
           obj.petFormat=Variable_initialize.readLine(gfileID,'PETFormat',commentSymbol,'string');
           obj.petPath=Variable_initialize.readLine(gfileID,'PETPath',commentSymbol,'string');
           obj.resultFormat=Variable_initialize.readLine(gfileID,'ResultFormat',commentSymbol,'string'); 
           obj.resultPath=Variable_initialize.readLine(gfileID,'ResultPath',commentSymbol,'string');
           obj.obsPath=Variable_initialize.readLine(gfileID,'OBSPath',commentSymbol,'string');
           
           
        
           
           
           
           
           
           obj.RainFactType=Variable_initialize.readLine(gfileID,'RainFactType',commentSymbol,'string');
           obj.RainFact=Variable_initialize.readLine(gfileID,'RainFact',commentSymbol,'double');
           obj.KsatType=Variable_initialize.readLine(gfileID,'KsatType',commentSymbol,'string');
           obj.Ksat=Variable_initialize.readLine(gfileID,'Ksat',commentSymbol,'double');
           obj.WMType=Variable_initialize.readLine(gfileID,'WMType',commentSymbol,'string');
           obj.WM=Variable_initialize.readLine(gfileID,'WM',commentSymbol,'double');
           obj.BType=Variable_initialize.readLine(gfileID,'BType',commentSymbol,'string');
           obj.B=Variable_initialize.readLine(gfileID,'B',commentSymbol,'double');
           obj.IMType=Variable_initialize.readLine(gfileID,'IMType',commentSymbol,'string');
           obj.IM=Variable_initialize.readLine(gfileID,'IM',commentSymbol,'double');
           obj.KEType=Variable_initialize.readLine(gfileID,'KEType',commentSymbol,'string');
           obj.KE=Variable_initialize.readLine(gfileID,'KE',commentSymbol,'double');
           obj.coeMType=Variable_initialize.readLine(gfileID,'coeMType',commentSymbol,'string');
           obj.coeM=Variable_initialize.readLine(gfileID,'coeM',commentSymbol,'double');
           
           %--------------read the landslide model input data-----------
           obj.Landslide_basicFormat=Variable_initialize.readLine(gfileID,'LandslideBasicFormat',commentSymbol,'string');
           obj.Landslide_basicPath=Variable_initialize.readLine(gfileID,'LandslideBasicPath',commentSymbol,'string');
           obj.UseMaxLocalNcores=Variable_initialize.readLine(gfileID,'UseMaxLocalNcores',commentSymbol,'boolean');
           obj.UserDefinedNcores=Variable_initialize.readLine(gfileID,'UserDefinedNcores',commentSymbol,'double');
           obj.Landslide_density=Variable_initialize.readLine(gfileID,'Landslide_density',commentSymbol,'double');
           obj.Divide_tile_number=Variable_initialize.readLine(gfileID,'Divide_tile_number',commentSymbol,'double');
           obj.min_ae=Variable_initialize.readLine(gfileID,'min_ae',commentSymbol,'double');
           obj.max_ae=Variable_initialize.readLine(gfileID,'max_ae',commentSymbol,'double');
           obj.min_be=Variable_initialize.readLine(gfileID,'min_be',commentSymbol,'double');
           obj.max_be=Variable_initialize.readLine(gfileID,'max_be',commentSymbol,'double');
           obj.LandslideOutStart_data = Variable_initialize.readLine(gfileID,'LandslideOutStart_data',commentSymbol,'string');
           obj.LandslideOutEnd_data = Variable_initialize.readLine(gfileID,'LandslideOutEnd_data',commentSymbol,'string');
           
           
           
           %--------------read the soil downscaling input data-----------
           obj.DownscalingBasicFormat = Variable_initialize.readLine(gfileID,'DownscalingBasicFormat',commentSymbol,'string');
           obj.DownscalingBasicPath = Variable_initialize.readLine(gfileID,'DownscalingBasicPath',commentSymbol,'string');
           
           
           %obj.expMType=Variable_initialize.readLine(gfileID,'expMType',commentSymbol,'string');
%            obj.coeRType=Variable_initialize.readLine(gfileID,'coeRType',commentSymbol,'string');
%            obj.coeSType=Variable_initialize.readLine(gfileID,'coeSType',commentSymbol,'string');
%            obj.KSType=Variable_initialize.readLine(gfileID,'KSType',commentSymbol,'string');
%            obj.KIType=Variable_initialize.readLine(gfileID,'KIType',commentSymbol,'string');
           
           

           %obj.expM=Variable_initialize.readLine(gfileID,'expM',commentSymbol,'double');
%            obj.coeR=Variable_initialize.readLine(gfileID,'coeR',commentSymbol,'double');
%            obj.coeS=Variable_initialize.readLine(gfileID,'coeS',commentSymbol,'double');
%            obj.KS=Variable_initialize.readLine(gfileID,'KS',commentSymbol,'double');
%            obj.KI=Variable_initialize.readLine(gfileID,'KI',commentSymbol,'double');
          
           
           
          
           
           obj.NOutPixs=Variable_initialize.readLine(gfileID,'NOutPixs',commentSymbol,'double');
           if obj.NOutPixs>0 
               obj.OutPixLats=zeros(obj.NOutPixs,1);
               obj.OutPixLons=obj.OutPixLats;
               obj.OutPixAreas=obj.OutPixLats;
               obj.OutPixNames=cell(obj.NOutPixs,1);
           end
           
           for i=1:obj.NOutPixs
               strname=Variable_initialize.readLine(gfileID,['OutPixName' num2str(i)],commentSymbol,'string');
               lat=Variable_initialize.readLine(gfileID,['OutPixLat' num2str(i)],commentSymbol,'double');
               lon=Variable_initialize.readLine(gfileID,['OutPixLon' num2str(i)],commentSymbol,'double');
               area=Variable_initialize.readLine(gfileID,['OutPixArea' num2str(i)],commentSymbol,'double');

               obj.OutPixLats(i)=lat;
               obj.OutPixLons(i)=lon;
               obj.OutPixAreas(i)=area;
               obj.OutPixNames{i}=strname;
           end
           
           obj.output_rain=Variable_initialize.readLine(gfileID,'GOVar_Rain',commentSymbol,'boolean');
           obj.output_epot=Variable_initialize.readLine(gfileID,'GOVar_EPot',commentSymbol,'boolean');
           obj.output_eact=Variable_initialize.readLine(gfileID,'GOVar_EAct',commentSymbol,'boolean');
           obj.output_W=Variable_initialize.readLine(gfileID,'GOVar_W',commentSymbol,'boolean');
           obj.output_SM=Variable_initialize.readLine(gfileID,'GOVar_SM',commentSymbol,'boolean');
           obj.output_runoff=Variable_initialize.readLine(gfileID,'GOVar_R',commentSymbol,'boolean');
           obj.output_rs=Variable_initialize.readLine(gfileID,'GOVar_RS',commentSymbol,'boolean');
           obj.output_ri=Variable_initialize.readLine(gfileID,'GOVar_RI',commentSymbol,'boolean');
%            obj.output_ExcS=Variable_initialize.readLine(gfileID,'GOVar_ExcS',commentSymbol,'boolean');
%            obj.output_ExcI=Variable_initialize.readLine(gfileID,'GOVar_ExcI',commentSymbol,'boolean');
           obj.output_infil=Variable_initialize.readLine(gfileID,'GOVar_Infil',commentSymbol,'boolean');
           obj.output_FS=Variable_initialize.readLine(gfileID,'GOVar_FS',commentSymbol,'boolean');
           obj.output_tot_infil=Variable_initialize.readLine(gfileID,'GOVar_tot_infil',commentSymbol,'boolean');
           obj.output_FS3D=Variable_initialize.readLine(gfileID,'GOVar_FS3D',commentSymbol,'boolean');
           obj.output_PF=Variable_initialize.readLine(gfileID,'GOVar_PF',commentSymbol,'boolean');
           obj.output_Volume=Variable_initialize.readLine(gfileID,'GOVar_Volume',commentSymbol,'boolean');
           obj.output_Area=Variable_initialize.readLine(gfileID,'GOVar_Area',commentSymbol,'boolean');
%            obj.output_pet=Variable_initialize.readLine(gfileID,'GOVar_PET',commentSymbol,'boolean');
           

        end  
    end

    methods(Static = true)
        function interval=CalTimeInterval(str,fmt)
            if strcmpi(fmt,'yyyy\DOY')
                content=strsplit(str,'\');
                interval=str2double(content{1})*365+str2double(content{2});
            else
                strOrg='';
                for i=1:length(fmt)
                    strOrg=strcat(strOrg,'0');
                end
                timeOrg=datenum(strOrg,fmt);
                time=datenum(str,fmt);
                interval=time-timeOrg;
            end
        end
        function value=readLine(gfileID,keyword,commentSymbol,type)
            value=-1;
            while value==-1
                tline = fgetl(gfileID);
                if strcmp(tline(1),commentSymbol)==1
                    continue;
                end
                strArr = regexp(tline,commentSymbol,'split');
                strContent=strArr{1};
                strContent=strtrim(strContent);
                if ~isempty(regexpi(strContent, keyword,'ONCE'))
                    strValue=regexp(strContent,'=','split');
                    switch type
                        case 'double'
                            value=str2double(strValue{2});
                        case 'string'
                            value=strrep(strValue{2},'"','');
                            value=strtrim(value);
                        case 'boolean'
                            value=Variable_initialize.yesno2boolean(strtrim(strValue{2}));
                    end
                else
                    %warning('warning: content disorder');
                end
            end
        end
        function bValue=yesno2boolean(value)
            if strcmpi(value,'yes')
               bValue=true;
           else
               bValue=false;
           end
        end
        function [dateForcStart,dateForcInter]=ForcingTimePar(strForcStart,strForcInter,fmt)
            if strcmpi(fmt,'yyyy\DOY')
                content=strsplit(strForcStart,'\');
                yearStart=str2double(content{1});
                doyStart=str2double(content{2});
                dateForcStart=datenum(yearStart,1,1);
                dateForcStart=dateForcStart+doyStart-1;
            else
                dateForcStart=datenum(strForcStart,fmt);
            end
            dateForcInter=Variable_initialize.CalTimeInterval(strForcInter,fmt);
        end
    end
end

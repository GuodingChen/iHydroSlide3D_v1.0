function [data, ref]=read_ascii_gis(filename)
            fid=fopen(filename,'r');
            headers=textscan(fid,'%s',6,'delimiter','\n');
            fclose(fid);

            for m=1:6
                strtmp=lower(char(headers{1}(m)));

                k=strfind(strtmp,'ncols');
                if ~isempty(k)
                   ref.ncols=str2double(strtmp(k+5:end));
                end

                k=strfind(strtmp,'nrows');
                if ~isempty(k)
                   ref.nrows=str2double(strtmp(k+5:end));
                end


                k=strfind(strtmp,'xllcorner');
                if ~isempty(k)
                   ref.xllcorner=str2double(strtmp(k+9:end));
                end

                k=strfind(strtmp,'yllcorner');
                if ~isempty(k)
                   ref.yllcorner=str2double(strtmp(k+9:end));
                end

                k=strfind(strtmp,'cellsize');
                if ~isempty(k)
                   ref.cellsize=str2double(strtmp(k+8:end));
                end

                k=strfind(strtmp,'nodata_value');
                if ~isempty(k)
                   ref.nodata=str2double(strtmp(k+12:end));
                end
            end
            fid=fopen(filename,'r');
            C=textscan(fid,'%f','headerlines',6);
            fclose(fid);
            tmp=C{1};
            data=reshape(tmp,[ref.ncols,ref.nrows])';
        end
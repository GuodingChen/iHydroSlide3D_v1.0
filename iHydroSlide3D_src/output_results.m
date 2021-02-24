function [success]=output_results(filename,format,geoinfo,projinfo,ref,data)
    success=0;
    switch lower(format)
        case {'asc','txt'}
            fid=fopen(filename,'w');
            fprintf(fid,'%-12s  %d\n','ncols', ref.ncols);
            fprintf(fid,'%-12s  %d\n','nrows', ref.nrows);
            fprintf(fid,'%-12s  %f\n','xllcorner', ref.xllcorner);
            fprintf(fid,'%-12s  %f\n','yllcorner', ref.yllcorner);
            fprintf(fid,'%-12s  %f\n','cellsize', ref.cellsize);
            fprintf(fid,'%-12s  %f\n','nodata', ref.nodata);

            for y=1:ref.nrows
                for x=1:ref.ncols
                    fprintf(fid,'%f ', data(y,x));
                end
                fprintf(fid,'\n');
            end
            fclose(fid);

        case 'tif'
            geotiffwrite(filename,single(data),geoinfo, 'GeoKeyDirectoryTag'...
            , projinfo.GeoTIFFTags.GeoKeyDirectoryTag...
            , 'TiffTags', struct('Compression','Deflate'));
        otherwise
            error(['The format of output files is ', format, ' that is not supported!']);
            success = 0;         
    end
    success=1;
end
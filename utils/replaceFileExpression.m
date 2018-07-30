function [] = replaceFileExpression(filepath, expression, replace)
    buffer = fileread(filepath);
    buffer = strrep(buffer, expression, replace);
    fid = fopen(filepath, 'w');
    fwrite(fid, buffer);
    fclose(fid);
end


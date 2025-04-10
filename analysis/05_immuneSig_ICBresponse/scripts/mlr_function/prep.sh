for file in ./*
do
    chmod 755 $file
    sed -i 's///g' $file
done


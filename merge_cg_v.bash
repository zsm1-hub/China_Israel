
# 设置默认值
Case="hf"
filt="uni"

# 根据Case设置文件名
if [ "$Case" = "smooth" ]; then
	    file='nowave'
    elif [ "$Case" = "hf" ]; then
	        file='wave'
	elif [ "$Case" = "HIT" ]; then
		    file='HIT2d'
	    else
		        echo "错误: 未知的 Case 值 '$Case'"
			    exit 1
fi

# 重命名文件
mv s2sflux_spec.0002.nc "s2sflux_spec_${Case}_${filt}.0002.nc"

# 合并文件
cdo merge \
	    "s2sflux_spec_all_${Case}.0002.nc" \
	        "../Parcels_data/${file}case_modified.nc" \
		    "../Parcels_data/${file}case_modified_cg_${filt}.nc"

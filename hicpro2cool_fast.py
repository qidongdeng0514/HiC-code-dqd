#!/usr/bin/env python3
import cooler
import numpy as np
import pandas as pd

def hicpro2cool(matrix_file, bed_file, output_cool):
    # 1. 读取BED文件，构建bins信息
    bins = pd.read_csv(bed_file, sep='\t', header=None, 
                       names=['chrom', 'start', 'end', 'bin_id'])
    
    # 2. 读取矩阵文件，构建像素表
    pixels = pd.read_csv(matrix_file, sep='\t', header=None,
                         names=['bin1_id', 'bin2_id', 'count'])
    # HiC-Pro矩阵是0-based索引，如果发现你的文件是1-based，减1即可
    # pixels[['bin1_id', 'bin2_id']] = pixels[['bin1_id', 'bin2_id']] - 1
    
    # 3. 创建Cooler文件
    cooler.create_cooler(
        cooler_uri=output_cool,
        bins=bins[['chrom', 'start', 'end']],
        pixels=pixels[['bin1_id', 'bin2_id', 'count']],
        dtypes={'count': np.float64}
    )
    print(f"✅ 已生成Cool文件: {output_cool}")

if __name__ == "__main__":
    import sys
    hicpro2cool(sys.argv[1], sys.argv[2], sys.argv[3])

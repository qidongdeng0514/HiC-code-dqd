#!/usr/bin/env python3
"""
hicpro2cool_fast.py - 将HiC-Pro输出的.matrix和.bed文件转换为.cool格式

使用方法:
    python hicpro2cool_fast.py <matrix_file> <bed_file> <output_cool>
"""

import cooler
import numpy as np
import pandas as pd
import sys
import os
import argparse
from pathlib import Path

def hicpro2cool(matrix_file, bed_file, output_cool, one_based=False):
    """
    将HiC-Pro的matrix和bed文件转换为cool格式
    """
    
    # 检查输入文件是否存在
    if not os.path.exists(matrix_file):
        raise FileNotFoundError(f"❌ 矩阵文件不存在: {matrix_file}")
    if not os.path.exists(bed_file):
        raise FileNotFoundError(f"❌ BED文件不存在: {bed_file}")
    
    print(f"📖 读取BED文件: {bed_file}")
    # 1. 读取BED文件，构建bins信息
    bins = pd.read_csv(bed_file, sep='\t', header=None, 
                       names=['chrom', 'start', 'end', 'bin_id'])
    print(f"   ✅ 共 {len(bins)} 个bins")
    
    print(f"📖 读取矩阵文件: {matrix_file}")
    # 2. 读取矩阵文件，构建像素表
    pixels = pd.read_csv(matrix_file, sep='\t', header=None,
                         names=['bin1_id', 'bin2_id', 'count'])
    
    # 如果matrix是1-based索引，转换为0-based
    if one_based:
        print("   🔄 将1-based索引转换为0-based...")
        pixels['bin1_id'] = pixels['bin1_id'] - 1
        pixels['bin2_id'] = pixels['bin2_id'] - 1
    
    # 过滤无效的bin索引
    max_bin = len(bins) - 1
    valid_mask = (pixels['bin1_id'] >= 0) & (pixels['bin1_id'] <= max_bin) & \
                 (pixels['bin2_id'] >= 0) & (pixels['bin2_id'] <= max_bin)
    if not valid_mask.all():
        invalid_count = (~valid_mask).sum()
        print(f"   ⚠️ 过滤掉 {invalid_count} 个无效像素")
        pixels = pixels[valid_mask]
    
    print(f"   ✅ 共 {len(pixels)} 个非零像素")
    print(f"   📊 最大count值: {pixels['count'].max():.0f}")
    
    # 确保数据类型正确
    pixels['bin1_id'] = pixels['bin1_id'].astype(np.int64)
    pixels['bin2_id'] = pixels['bin2_id'].astype(np.int64)
    pixels['count'] = pixels['count'].astype(np.float64)
    
    # 3. 创建Cooler文件
    print(f"💾 创建cool文件: {output_cool}")
    output_cool = Path(output_cool)
    output_cool.parent.mkdir(parents=True, exist_ok=True)
    
    # 计算分辨率，确保是Python原生int类型
    resolution = int(bins['end'].iloc[0] - bins['start'].iloc[0])
    
    # 使用位置参数，避免JSON序列化问题
    cooler.create_cooler(
        str(output_cool),  # 第一个参数是URI，位置参数
        bins=bins[['chrom', 'start', 'end']],
        pixels=pixels[['bin1_id', 'bin2_id', 'count']],
        dtypes={'count': np.float64},
        assembly='hg38'
    )
    
    print(f"✅ 成功生成Cool文件: {output_cool}")
    
    # 验证生成的cool文件
    try:
        clr = cooler.Cooler(str(output_cool))
        print(f"   📊 分辨率: {clr.binsize} bp")
        print(f"   📊 Bin数量: {clr.info['nbins']:,}")
        print(f"   📊 像素数量: {clr.info['nnz']:,}")
    except Exception as e:
        print(f"   ⚠️ 验证失败: {e}")

def main():
    parser = argparse.ArgumentParser(
        description="将HiC-Pro输出的.matrix和.bed文件转换为.cool格式",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:
  %(prog)s A375_ActD_merged_5000.matrix A375_ActD_merged_5000_abs.bed A375_ActD_merged_5000.cool
  %(prog)s --one-based A375_ActD_merged_5000.matrix A375_ActD_merged_5000_abs.bed A375_ActD_merged_5000.cool
        """
    )
    
    parser.add_argument(
        'matrix_file',
        type=str,
        help='HiC-Pro输出的.matrix文件路径'
    )
    parser.add_argument(
        'bed_file',
        type=str,
        help='HiC-Pro输出的_abs.bed文件路径'
    )
    parser.add_argument(
        'output_cool',
        type=str,
        help='输出的.cool文件路径'
    )
    parser.add_argument(
        '--one-based', '-1',
        action='store_true',
        help='如果.matrix文件是1-based索引（默认是0-based）'
    )
    parser.add_argument(
        '--version', '-v',
        action='version',
        version='hicpro2cool_fast.py 1.0.2'
    )
    
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)
    
    args = parser.parse_args()
    
    try:
        hicpro2cool(
            matrix_file=args.matrix_file,
            bed_file=args.bed_file,
            output_cool=args.output_cool,
            one_based=args.one_based
        )
    except Exception as e:
        print(f"❌ 转换失败: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()

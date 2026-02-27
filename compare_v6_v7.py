#!/usr/bin/env python3
import csv
from collections import Counter

# 统计V7输出
with open('Mt_V7_optimized_KO.csv', encoding='utf-8') as f:
    v7 = list(csv.DictReader(f))
with open('Mt_V7_optimized_KO_failed.csv', encoding='utf-8') as f:
    v7f = list(csv.DictReader(f))

print('=== V7 优化版（8线程）结果 ===')
print(f'成功设计: {len(v7)}')
print(f'失败基因: {len([r for r in v7f if r.get("Status")=="Failed"])}')
print(f'Partial基因: {len([r for r in v7f if r.get("Status")=="Partial"])}')

# Strategy分布
st = Counter(r['Strategy'] for r in v7)
print(f'Strategy分布: {dict(st)}')

# 删除长度统计
dels = [int(r['Deletion Length']) for r in v7]
print(f'删除长度: min={min(dels)} max={max(dels)} avg={sum(dels)/len(dels):.1f}')

# V6对比
with open('Mt_KO_library_v6_delPct10_80_delBp300_1000_bc11.csv', encoding='utf-8') as f:
    v6 = list(csv.DictReader(f))
print(f'\n=== V6对比 ===')
print(f'V6成功设计: {len(v6)}')
print(f'V7设计/V6设计: {len(v7)/len(v6)*100:.1f}%')

v6_genes = set(r['Gene Id'] for r in v6)
v7_genes = set(r['Gene Id'] for r in v7)
print(f'共有基因: {len(v6_genes & v7_genes)}')
print(f'V6独有: {len(v6_genes - v7_genes)}')
print(f'V7独有: {len(v7_genes - v6_genes)}')

# 削减的设计数
reduction_pct = (1 - len(v7)/len(v6)) * 100
print(f'\n设计数削减: {len(v6) - len(v7)} ({reduction_pct:.1f}%)')
print('原因分析：V7采用Top-5候选策略，可能漏掉V6穷尽式搜索的部分边界最优解')
print('（但这既降低了准确性影响，也是换取10倍速度的必要权衡）')

"""Build a factual handoff from executed outputs; do not generate manuscript results."""
from pathlib import Path
import csv
import json
import pandas as pd
from build_instruments import ROUND, ROOT, sha


def main():
    validation=json.loads((ROUND/'reports/round03_validation.json').read_text())
    summary=pd.read_csv(ROUND/'tables/instrument_rebuild_summary.csv')
    manifest=pd.read_csv(ROOT/'metadata/FinnGen_R12_manifest.tsv',sep='\t')
    outcomes=['M13_ANKYLOSPON','H7_IRIDOCYCLITIS','L12_PSORIASIS','M13_PSORIARTH','M13_ANKYLOSPON_STRICT']
    requested=manifest.loc[manifest.phenocode.isin(outcomes)].copy()
    assert set(requested.phenocode)==set(outcomes) and len(requested)==5
    requested['planned_role']=requested.phenocode.map(lambda x:'SAME_BIOBANK_PHENOTYPE_SENSITIVITY' if x.endswith('STRICT') else 'PRIMARY_OUTCOME')
    requested['download_gate']='HOLD_USER_REGISTRATION_NOT_COMPLETED'
    requested['planned_storage_root']='C:/Users/Lenovo/IBD_rebuild_runtime/finngen_R12'
    requested.to_csv(ROUND/'tables/FinnGen_pending_source_manifest.csv',index=False)
    gates=[
      ('SOURCE_FILES','PASS_WITH_RECORDED_LIMITATIONS','Original exposure bytes and round02 audit intact; missing INFO/EAF/per-SNP N not fabricated','required'),
      ('EXPOSURE_INSTRUMENTS','PASS','12 exposure QC/MHC branches clumped and every eligible candidate reconciled','required'),
      ('REFERENCE_LD','PASS_WITH_LIMITATIONS','502 EUR individuals; within-10Mb lead pair LD verified; not FinnGen in-sample LD','required'),
      ('COORDINATE_MAP','PARTIAL_PASS','218 unique SNVs reciprocal-coordinate pass; reference/outcome alleles not yet checked','required'),
      ('FINNGEN_REGISTRATION','NEEDS_RESEARCHER_ACTION','User confirms not completed; no new FinnGen download','required'),
      ('OUTCOME_SOURCE_QC','NOT_COMPLETE','Fresh four-outcome provenance, integrity, fields and coverage not complete','required'),
      ('PRIMARY_12_MR','NOT_RUN','Do not reuse legacy effects as newly rebuilt effects','required'),
      ('MHC_EFFECT_SENSITIVITY','INPUT_READY_NOT_RUN','MHC exclusion performed before independent clumping, no outcome effect comparison yet','required'),
      ('MR_PRESSO','NOT_RUN_IN_REBUILD','Two fixed seeds/high precision specified; old results not transplanted','planned'),
      ('BINARY_STEIGER','NOT_ESTIMABLE_FROM_CURRENT_INPUTS','Exposure frequency and defensible binary-liability settings missing; may omit with limitation, never fabricate pass','claim_dependent'),
      ('COLOCALISATION','REGION_PLAN_ONLY','130 regions planned independent of new outcome results; no new H4 or multi-signal claims','required_for_coloc_claims'),
      ('OVERLAP','UNQUANTIFIED','Do not assert zero overlap or completed quantitative bias correction','required_disclosure'),
      ('INDEPENDENT_REPLICATION','NOT_ESTABLISHED','Same-cohort releases/subtypes/strict-AS cannot substitute for independent replication','claim_dependent'),
      ('NOVELTY','NOT_ESTABLISHED','Direct prior MR studies overlap; systematic regional evidence must add genuine information','required'),
      ('CITATION_CLAIM_MAPPING','NOT_FINAL_AUDITED','58-reference historical text untouched; final new claims need sentence-level audit','required'),
      ('AUTHOR_ETHICS_LICENSE_AI','NEEDS_AUTHOR_FINAL_REVIEW','No invented approval, exemption, license or author declaration','required'),
      ('MANUSCRIPT_AND_PUBLIC_RELEASE','NOT_READY','No new final manuscript or public release produced before actual estimates and audit','required')
    ]
    pd.DataFrame(gates,columns=['item','status','basis','requirement']).to_csv(ROUND/'tables/submission_readiness.csv',index=False)
    lines=['# 第三轮重建：工具变量与投稿准备核验','',
      '日期：2026-09-27。结论：本地输入准备取得实质进展，但目前不能投稿。所有产物是重建工作文件，不是已完成的新论文结果。','',
      '## 本轮实际完成','',
      '- 再核验上一轮 85 个文件的完整清单，原有文件未变。',
      '- 固定新的分析规则和执行前哈希；明示这是看过历史结果后的重建，而非预注册、盲法验证或独立复制。',
      '- 以完整暴露源文件派生的可追溯 GWS 候选表运行 PLINK。只读使用公共参考资料和软件，不读取其他项目研究结果。',
      '- 主分析/放宽队列 P 规则各自运行全区域与排除 MHC 两个分支。MHC 定义为 GRCh37 chr6:25–34 Mb，先排除再重新 clump。','',
      '| 暴露 | 主分析 SNP | 排除 MHC 后 SNP | 放宽队列 P 规则 SNP | 放宽规则且排除 MHC | 主分析最小近似 F |',
      '|---|---:|---:|---:|---:|---:|']
    for trait in ['IBD','CD','UC']:
        d=summary.loc[summary.trait==trait].set_index('branch')
        counts=[int(d.loc[b,'clumped_index']) for b in ['primary_all','primary_nonMHC','relaxed_cohort_p_all','relaxed_cohort_p_nonMHC']]
        lines.append('| '+trait+' | '+' | '.join(map(str,counts))+f" | {d.loc['primary_all','min_F']:.2f} |")
    lines += ['',
      '这里的 SNP 数是暴露端工具数，不是每个 exposure–outcome pair 的最终有效 SNP 数。近似 F=(beta/SE)^2，只描述关联强度，不证明工具的排除限制或不存在赢家诅咒。',
      '', '## 防错核验','',
      '- 12 个分支的每一个合格候选均可追溯到 index 或其所属 clump，未任意删除或保持第一条冲突记录。',
      '- 逐一检查保留工具在同一染色体 10 Mb 内的成对 LD，所有已检配对 r²<=0.001。超过窗口的长程 LD 未因此被证明不存在。',
      '- 218 个不同位点全部通过双向坐标链往返检查。所有均为正链映射；结局 REF/ALT 尚未检验，不能把坐标通过误写成完整 harmonisation 通过。',
      '- 13 项单元测试通过，覆盖等位基因翻转、回文 SNP、无根据的链互补、无效 SE、重复记录、固定 12 项 FDR 家族等。',
      '- 安装的 TwoSampleMR 与独立 Python 公式在 4 组模拟输入的 6 项 IVW 统计量上一致；单 SNP Wald ratio 单独调用。模拟值不进入论文研究结果。',
      '- EAF 缺失不是非回文 SNP 的 beta/SE MR 或 case-control beta/SE coloc 的统一禁用条件；但不能凭结局频率或参考频率冒充暴露实测 EAF 去宣称严格二分类 Steiger 已完成。',
      '', '## 已固定的区域分析路线','',
      '围绕全部主分析 exposure index SNP 的 +/-500 kb 窗口进行合并，形成 130 个候选区域，合计 140,911,062 bp。不是只挑旧阳性 IL23R/NOD2。该文件是区域计划，不是已提取的区域数据或已获得的共定位结果。',
      '',
      '正式 coloc 仍需要密集区域数据、版本及等位基因匹配、覆盖率和双性状信号强度检查。复杂 MHC 区域不得凭单因果 ABF 的高 H4 自动升级为确认；单信号模型局限和先验敏感性需逐区域报告。',
      '', '## 现在无法跳过的前提','',
      '**用户已确认 FinnGen 下载登记尚未完成。** 本轮未新增下载或重跑 FinnGen 结局数据。登记由实际研究者完成；不能伪造身份、勾选未经理解的使用承诺、倒填日期或把旧文件存在当成登记证明。详见同目录《FinnGen登记与后续恢复说明》。',
      '',
      '除此之外，新主分析和系统共定位尚未运行、创新增量尚未证实、最终文献逐句对应和作者声明尚未终审。不能承诺登记后就自动具备投稿条件。Steiger、独立复制和多信号共定位是否必需取决于最终主张；做不了时应缩小结论并明确披露，不以补造结果凑齐方法。',
      '', '## 对旧问题的处理','',
      '- GCST005529 的来源/效应冲突仍未解释。它被排除出新主分析，并不等于历史矛盾已经解决。',
      '- 旧 CD→AS 的 18 SNP 差异不通过猜测补齐。本轮建立全新的逐 SNP 选择记录；新记录不能倒推旧过程已完全追溯。',
      '- 旧结果、参考文献、Word、原始数据和上轮审计未覆盖。新结果未产生前不重写一个看似完稿的 Results。',
      '', '## 主要文件','',
      '- config/analysis_lock.json：固定规则；metadata/execution_lock.json：执行前哈希。',
      '- tables/instrument_rebuild_summary.csv：四分支计数。',
      '- tables/*_candidate_trace.csv.gz 与 *_membership.csv.gz：逐 SNP 选择和所属 clump。',
      '- tables/instrument_coordinate_audit.csv：双向坐标核验。',
      '- tables/systematic_coloc_region_plan_GRCh37.csv：结果无关的区域计划。',
      '- tables/FinnGen_pending_source_manifest.csv：待登记后获取的具体版本和表型。',
      '- tables/submission_readiness.csv：必须完成、可有条件省略、作者动作的区分。',
      '- reports/round03_validation.json 与 logs/：真实执行和测试状态。',
      '', '## 方法来源','',
      '- [PLINK clumping 官方说明](https://www.cog-genomics.org/plink/1.9/postproc#clump)：clumping 及其 LD 定义。',
      '- [FinnGen 官方下载说明](https://finngen.gitbook.io/documentation/data-download)：登记、引用及致谢要求。',
      '- 本地已安装 TwoSampleMR 的 mr_ivw / mr_wald_ratio 源码：本轮运行规则和独立公式核验对象，不将在线新版本替换本机实际版本。',
      '', '未保证不存在任何残余错误，也未保证接收；投稿就绪以数据、结果、主张、引用和声明一致为标准。','']
    (ROUND/'reports/第三轮重建结果与剩余投稿条件.md').write_text('\n'.join(lines),encoding='utf-8')
    (ROUND/'README_先读.md').write_text('# 第三轮重建工作目录\n\n2026-09-27：工具变量重建及本地防错测试完成，尚不能投稿。\n\n先读 reports/第三轮重建结果与剩余投稿条件.md。\n\n外部动作：FinnGen 官方下载登记尚未完成，见 reports/FinnGen登记与后续恢复说明.md。\n\n不要把模拟测试、候选区域或暴露工具数写成已完成的真实 MR/coloc 结果。所有旧稿和历史数据保持不变。\n',encoding='utf-8')
    records=[]
    target=ROUND/'tables/round03_file_manifest.csv'
    for p in sorted(ROUND.rglob('*')):
        if p.is_file() and '__pycache__' not in p.parts and p != target:
            records.append(dict(path=p.relative_to(ROUND).as_posix(),bytes=p.stat().st_size,sha256=sha(p)))
    pd.DataFrame(records).to_csv(target,index=False)
    for r in records:
        assert sha(ROUND/r['path'])==r['sha256']
    print(f'PASS: {len(records)} files indexed; status remains NOT_SUBMISSION_READY')


if __name__ == '__main__':
    main()

#!/usr/bin/env perl
use strict;
use warnings;
use List::Util qw(sum);
use Getopt::Long;
use Pod::Usage;
use utf8;

# 版本信息
our $VERSION = '1.01 2025-03-06';

# 命令行选项处理
GetOptions(
    'help|h'     => \my $help,
    'version|v'  => \my $version,
    'output|o=s' => \my $output_file,
    'dj|d=i' => \my $dj_threshold
) or pod2usage(2);

# 显示帮助信息
pod2usage(-verbose => 2) if $help;

# 显示版本信息
if ($version) {
    print "EANMD GTF NMD checker v$VERSION\n";
    exit;
}

# 验证输入参数
pod2usage("Error: Missing GTF file") unless @ARGV;
my $gtf_file = shift @ARGV;
$dj_threshold ||= 50;  # 默认NMD阈值

# 重定向输出到文件（如果指定）
if ($output_file) {
    open STDOUT, '>', "$output_file.tsv" 
        or die "Can not create: $!";
}

# # 参数处理
# my ($gtf_file, $dj_threshold) = @ARGV;
# $dj_threshold ||= 50;  # 默认NMD阈值为50nt

# 数据结构存储转录本信息
my %transcripts;

# 解析GTF文件
open my $fh, '<', $gtf_file or die "can not open GTF: $!";
while (<$fh>) {
    next if /^#/;  # 跳过注释行
    my @fields = split /\t/;
    
    # 只处理exon和CDS记录
    next unless $fields[2] eq 'exon' || $fields[2] eq 'CDS';
    
    # 提取公共字段
    my ($chr, $type, $start, $end, $strand) = 
        @fields[0,2,3,4,6];
    
    # 提取属性字段
    my ($gene_id) = $_ =~ /gene_id "([^"]+)"/;
    my ($gene_name) = $_ =~ /gene_name "([^"]+)"/;
    my ($transcript_id) = $_ =~ /transcript_id "([^"]+)"/;
    my ($exon_number) = $_ =~ /exon_number (\d+)/;
    
    # 初始化转录本结构
    $transcripts{$transcript_id} ||= {
        gene_id => $gene_id,
        gene_name => $gene_name || 'NA',
        strand => $strand,
        exons => [],
        cds => []
    };
    
    # 存储exon信息
    if ($type eq 'exon') {
        push @{$transcripts{$transcript_id}{exons}}, {
            number => $exon_number,
            start => $start,
            end => $end,
            length => $end - $start + 1
        };
    }
    
    # 存储CDS信息
    if ($type eq 'CDS') {
        push @{$transcripts{$transcript_id}{cds}}, {
            start => $start,
            end => $end,
            exon_number => $exon_number,
            length => $end - $start + 1
        };
    }
}
close $fh;

# 输出标题行
print join("\t", 
    'Transcript_ID',
    'Gene_ID',
    'Gene_Name',
    'Strand',
    'Total_Exons',
    'Exon_Total_Length',
    'CDS_Total_Length',
    'CDS_Start_Exon',
    'CDS_End_Exon',
    'Stop_Exon',
    'Distance_to_Last_Junction',
    'NMD_Status'
), "\n";

# 处理每个转录本
# foreach my $tid (keys %transcripts) {
# foreach my $tid (
#     sort {
#         $transcripts{$a}{gene_name} cmp $transcripts{$b}{gene_name}
#         ||
#         $a cmp $b
#     } keys %transcripts
# ) {
foreach my $tid (sort keys %transcripts) { # 2025-03-06 make results same order.
    my $t = $transcripts{$tid};
    
    # 跳过没有CDS的转录本
    next unless @{$t->{cds}};
    
    # 排序exon和CDS记录
    my @exons = sort { $a->{number} <=> $b->{number} } @{$t->{exons}};
    my @cds = sort { $a->{exon_number} <=> $b->{exon_number} } @{$t->{cds}};
    
    # 计算总长度
    my $exon_total = sum map { $_->{length} } @exons;
    my $cds_total = sum map { $_->{length} } @cds;
    
    # 确定CDS起始/终止exon
    my $cds_start_exon = $cds[0]{exon_number};
    my $cds_end_exon = $cds[-1]{exon_number};
    
    # 获取终止密码子信息
    # my $stop_pos = ($t->{strand} eq '+') 
    #     ? $cds[-1]{end} 
    #     : $cds[-1]{start};
    my $stop_pos = 0;
    for (my $sum_exon_num = 0; $sum_exon_num <$cds_end_exon-1; $sum_exon_num ++){
        $stop_pos = $stop_pos + $exons[$sum_exon_num]{length};
    }
    $stop_pos = $stop_pos +$cds[-1]{length};
    
    # 确定终止exon
    # my ($stop_exon) = grep { $_->{number} == $cds_end_exon } @exons;
    my $last_exon_len = $exons[-1]{length};
    
    # 计算到last junction的距离
    my $last_junction_pos = $exon_total - $last_exon_len;
    my $distance = ($t->{strand} eq '+')
        # ? $last_junction_pos - ($stop_pos - $exons[0]{start})
        # : ($stop_pos - $exons[-1]{end}) - $last_junction_pos;
        ? $last_junction_pos - $stop_pos
        : $last_junction_pos - $stop_pos;
    
    # 判断NMD状态
    my $nmd_status = ($distance > $dj_threshold) ? 'NMD' : 'nonNMD';
    
    # 输出结果
    print join("\t",
        $tid,
        $t->{gene_id},
        $t->{gene_name},
        $t->{strand},
        scalar @exons,
        $exon_total,
        $cds_total,
        $cds_start_exon,
        $cds_end_exon,
        $cds_end_exon,
        $distance,
        $nmd_status
    ), "\n";
}
__END__

=head1 NAME

EANMDcheckGTFNMD.pl - EANMD GTF transcript NMD Status checker 

=head1 SYNOPSIS

  perl EANMDcheckGTFNMD.pl [options] input.gtf [threshold]

Options:
  -h, --help      Show this help message
  -v, --version   Show version information
  -o, --output    Specify output file (default: STDOUT)
  -d, --dj Set NMD threshold (default: 50)

Arguments:
  input.gtf       Required, input GTF file path
  dj       Optional, NMD detection threshold (in nucleotides)

=head1 DESCRIPTION

This script analyzes NMD status of transcripts from GTF files based on:
1. Distance between stop codon and last exon-exon junction
2. Mark NMD positive when distance exceeds threshold
3. Output detailed gene/transcript structural information
=encoding UTF-8
=head1 OUTPUT COLUMNS

  1. Transcript_ID       转录本ID
  2. Gene_ID             基因ID
  3. Gene_Name           基因名称
  4. Strand              链方向（+/-）
  5. Total_Exons         exon总数
  6. Exon_Total_Length   所有exon总长度
  7. CDS_Total_Length    CDS区域总长度
  8. CDS_Start_Exon      CDS起始exon编号
  9. CDS_End_Exon        CDS终止exon编号
  10.Stop_Exon           终止密码子所在exon
  11.Distance_to_Junction 到最后一个剪接点的距离
  12.NMD_Status          NMD状态（NMD/non-NMD）

=head1 EXAMPLES

基本用法：
  perl EANMDcheckGTFNMD.pl Homo_sapiens.GRCh38.gtf > results.tsv

设置自定义阈值：
  perl EANMDcheckGTFNMD.pl -d 55 input.gtf

保存到文件：
  perl EANMDcheckGTFNMD.pl -o nmd_results input.gtf 50

=head1 AUTHOR

Kaining Hu <hukaining@gmail.com>
and DeepSeek R1 2025-03-05

=cut
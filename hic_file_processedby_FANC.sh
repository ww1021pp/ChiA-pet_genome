######this scipts is to processed the .hic file produced by CLoops and juicer software to the .hic file can be used in FAN-C
#(https://vaquerizaslab.github.io/fanc/fanc-executable/fanc_basic.html)###
/usr/bin/bash
for i in $(ls S*_hic1.hic)
do
echo $i
sample=`basename $i _hic1.hic`
fanc hic -m ICE -n ${i}@100kb ./fanc_hic/${sample}_iced_fromjuicer_100kb.hic -b 100000 --deepcopy
fanc hic -m ICE -n ${i}@50kb ./fanc_hic/${sample}_iced_fromjuicer_50kb.hic -b 50000 --deepcopy
fanc hic -m ICE -n ${i}@5kb ./fanc_hic/${sample}_iced_fromjuicer_5kb.hic -b 5000 --deepcopy
done


fanc directionality S1_hic1.hic@10kb ./fanc_domain/S1_fanc_10kb.directionality -w 100000 150000 200000 2000000 3000000
fancplot -o ./fanc_domain/fanc_example_10kb_tads_directionality_juicer.png 7:2mb-10mb \
     -p triangular -c Reds S1_hic1.hic@10kb -m 2mb -vmin 0 -vmax 0.05 \
     -p scores ./fanc_domain/S1_fanc_10kb.directionality
     
     
 fanc hic -m ICE S1_hic1.hic@10kb S1_juiced_iced_norm_10kb.hic -b 10kb --deepcopy  ##normlized juicer .hic file with a special resolution
 fanc directionality S1_juiced_iced_norm_10kb.hic ./fanc_domain/S1_10kb_juicerNorm.directionality -w 100000 150000 200000 2000000 3000000
      fancplot -o ./fanc_domain/fanc_example_10kb_tads_directionality_juicer_norm.png 7:2mb-10mb \
     -p triangular -c Reds S1_juiced_iced_norm_10kb.hic -m 2mb -vmin 0 -vmax 0.05 \
     -p scores ./fanc_domain/S1_10kb_juicerNorm.directionality
  
  fanc directionality ./fanc_hic/S1_fanc_10kb.hic ./fanc_domain/S1_10kb_fancraw.directionality -w 100000 150000 200000 2000000 3000000
     fancplot -o ./fanc_domain/fanc_example_10kb_tads_directionality_fanc_raw.png 7:2mb-10mb \
     -p triangular -c Reds ./fanc_hic/S1_fanc_10kb.hic -m 2mb -vmin 0 -vmax 0.05 \
     -p scores ./fanc_domain/S1_10kb_fancraw.directionality
     
     
     
   fanc directionality ./fanc_hic/S1_fanc_iced_norm_10kb.hic ./fanc_domain/S1_10kb_fancNorm.directionality -w 100000 150000 200000 2000000 3000000
     
     fancplot -o ./fanc_domain/fanc_example_10kb_tads_directionality_fancNorm.png 7:2mb-10mb \
     -p triangular -c Reds ./fanc_hic/S1_fanc_iced_norm_10kb.hic -m 2mb -vmin 0 -vmax 0.05 \
     -p scores ./fanc_domain/S1_10kb_fancNorm.directionality    
     
     

#include "pfp_meta.hpp"
#include <memory>
#include <cstdlib>
namespace {pm::Meta meta;std::unique_ptr<pm::Map>keep;}
extern "C" void sxgc_joint_filter_begin(void){try{
    auto path=getenv("SXGC_JOINT_KEEP"),prefix=getenv("SXGC_JOINT_MARK_PREFIX");
    if(!path&&!prefix)return;if(!path||!prefix)pm::fail("joint filter requires both keep bitmap and mark metadata");
    meta=pm::Meta::read(std::string(prefix)+".meta");keep=std::make_unique<pm::Map>(path);
    if(keep->size!=(meta.n+7)/8)pm::fail("joint keep bitmap dimensions mismatch");
}catch(const std::exception&e){fprintf(stderr,"joint filter failed: %s\n",e.what());exit(1);}}
extern "C" int sxgc_joint_filter_keep(int contig,long long anchor,int reverse){
    if(!keep)return 1;if(contig<0||size_t(contig)>=meta.contigs.size()||anchor<0){fputs("joint filter coordinate out of range\n",stderr);exit(1);}
    auto&c=meta.contigs[contig];uint64_t p=anchor;
    if(p>c.length||c.length-p<12||c.forward==UINT64_MAX){fputs("joint filter invalid sampled anchor\n",stderr);exit(1);}
    uint64_t at=reverse?c.reverse+c.length-p-12:c.forward+p;
    return keep->data[at>>3]&(1u<<(at&7));
}

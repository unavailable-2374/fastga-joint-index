// Included inside sg. Each stream owns its cache; clones never share mutable
// state across worker threads. Cursors pin immutable blocks while using them.
template<class T,size_t Capacity=4> class ReaderBlockCache {
    struct Slot { uint64_t block=UINT64_MAX,stamp=0; std::shared_ptr<const T> value; };
    std::array<Slot,Capacity> slots_{};
    uint64_t clock_=0;
public:
    uint64_t hits=0,misses=0;
    template<class Load> std::shared_ptr<const T> get(uint64_t block,Load load) {
        for(auto &slot:slots_)if(slot.block==block && slot.value){
            slot.stamp=++clock_;++hits;return slot.value;
        }
        auto value=std::make_shared<const T>(load()); // Publish only validated data.
        auto victim=std::min_element(slots_.begin(),slots_.end(),[](const Slot&a,const Slot&b){return a.stamp<b.stamp;});
        *victim={block,++clock_,value};++misses;return value;
    }
};

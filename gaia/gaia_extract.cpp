/**
 * @file
 * @brief Extract Data from GAIR EDR3 gaia_source/*.gz.csv
 * @author AKIYAMA Kouhei
 */
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <string>
#include <stdexcept>
#include <vector>
#include <utility>
#include <limits>
#include <filesystem>
#include <optional>
#include <zlib.h>

// gaia_source csv column indices:
const int col_solution_id = 0;
const int col_designation = 1;
const int col_source_id = 2;
const int col_random_index = 3;
const int col_ref_epoch = 4;
const int col_ra = 5;
const int col_ra_error = 6;
const int col_dec = 7;
const int col_dec_error = 8;
const int col_parallax = 9;
const int col_parallax_error = 10;
const int col_parallax_over_error = 11;
const int col_pm = 12;
const int col_pmra = 13;
const int col_pmra_error = 14;
const int col_pmdec = 15;
const int col_pmdec_error = 16;
const int col_ra_dec_corr = 17;
const int col_ra_parallax_corr = 18;
const int col_ra_pmra_corr = 19;
const int col_ra_pmdec_corr = 20;
const int col_dec_parallax_corr = 21;
const int col_dec_pmra_corr = 22;
const int col_dec_pmdec_corr = 23;
const int col_parallax_pmra_corr = 24;
const int col_parallax_pmdec_corr = 25;
const int col_pmra_pmdec_corr = 26;
const int col_astrometric_n_obs_al = 27;
const int col_astrometric_n_obs_ac = 28;
const int col_astrometric_n_good_obs_al = 29;
const int col_astrometric_n_bad_obs_al = 30;
const int col_astrometric_gof_al = 31;
const int col_astrometric_chi2_al = 32;
const int col_astrometric_excess_noise = 33;
const int col_astrometric_excess_noise_sig = 34;
const int col_astrometric_params_solved = 35;
const int col_astrometric_primary_flag = 36;
const int col_nu_eff_used_in_astrometry = 37;
const int col_pseudocolour = 38;
const int col_pseudocolour_error = 39;
const int col_ra_pseudocolour_corr = 40;
const int col_dec_pseudocolour_corr = 41;
const int col_parallax_pseudocolour_corr = 42;
const int col_pmra_pseudocolour_corr = 43;
const int col_pmdec_pseudocolour_corr = 44;
const int col_astrometric_matched_transits = 45;
const int col_visibility_periods_used = 46;
const int col_astrometric_sigma5d_max = 47;
const int col_matched_transits = 48;
const int col_new_matched_transits = 49;
const int col_matched_transits_removed = 50;
const int col_ipd_gof_harmonic_amplitude = 51;
const int col_ipd_gof_harmonic_phase = 52;
const int col_ipd_frac_multi_peak = 53;
const int col_ipd_frac_odd_win = 54;
const int col_ruwe = 55;
const int col_scan_direction_strength_k1 = 56;
const int col_scan_direction_strength_k2 = 57;
const int col_scan_direction_strength_k3 = 58;
const int col_scan_direction_strength_k4 = 59;
const int col_scan_direction_mean_k1 = 60;
const int col_scan_direction_mean_k2 = 61;
const int col_scan_direction_mean_k3 = 62;
const int col_scan_direction_mean_k4 = 63;
const int col_duplicated_source = 64;
const int col_phot_g_n_obs = 65;
const int col_phot_g_mean_flux = 66;
const int col_phot_g_mean_flux_error = 67;
const int col_phot_g_mean_flux_over_error = 68;
const int col_phot_g_mean_mag = 69;
const int col_phot_bp_n_obs = 70;
const int col_phot_bp_mean_flux = 71;
const int col_phot_bp_mean_flux_error = 72;
const int col_phot_bp_mean_flux_over_error = 73;
const int col_phot_bp_mean_mag = 74;
const int col_phot_rp_n_obs = 75;
const int col_phot_rp_mean_flux = 76;
const int col_phot_rp_mean_flux_error = 77;
const int col_phot_rp_mean_flux_over_error = 78;
const int col_phot_rp_mean_mag = 79;
const int col_phot_bp_n_contaminated_transits = 80;
const int col_phot_bp_n_blended_transits = 81;
const int col_phot_rp_n_contaminated_transits = 82;
const int col_phot_rp_n_blended_transits = 83;
const int col_phot_proc_mode = 84;
const int col_phot_bp_rp_excess_factor = 85;
const int col_bp_rp = 86;
const int col_bp_g = 87;
const int col_g_rp = 88;
const int col_dr2_radial_velocity = 89;
const int col_dr2_radial_velocity_error = 90;
const int col_dr2_rv_nb_transits = 91;
const int col_dr2_rv_template_teff = 92;
const int col_dr2_rv_template_logg = 93;
const int col_dr2_rv_template_fe_h = 94;
const int col_l = 95;
const int col_b = 96;
const int col_ecl_lon = 97;
const int col_ecl_lat = 98;


unsigned int numErrors = 0;
template<typename... Args>
void error(const char *fmt, Args... args)
{
    std::fprintf(stderr, fmt, args...);
    std::fprintf(stderr, "\n");
    ++numErrors;
}


class FloatOutputFile
{
    std::unique_ptr<float[]> buffer_;
    const std::size_t bufferCapacity_;
    std::size_t bufferSize_;
    std::FILE *file_;
public:
    FloatOutputFile(const char *filename, std::size_t bufferCapacity)
        : buffer_(new float[bufferCapacity])
        , bufferCapacity_(bufferCapacity)
        , bufferSize_(0)
        , file_(std::fopen(filename, "wb"))
    {}
    ~FloatOutputFile()
    {
        if(file_){
            flush();
            std::fclose(file_);
        }
    }
    bool isOpen() const {return file_ != NULL;}
    void flush()
    {
        if(bufferSize_ > 0){
            const std::size_t blocksWritten = std::fwrite(buffer_.get(), bufferSize_ * sizeof(float), 1, file_);
            if(blocksWritten != 1){
                error("fwrite error");
            }
            bufferSize_ = 0;
        }
    }
    void write(float value)
    {
        if(bufferSize_ >= bufferCapacity_){
            flush();
        }
        buffer_[bufferSize_++] = value;
    }
};

class ByteBuffer
{
    std::unique_ptr<char[]> bytes_;
    std::size_t capacity_;
    std::size_t size_;
public:
    ByteBuffer()
        : bytes_()
        , capacity_(0)
        , size_(0)
    {}
    char *get() {return bytes_.get();}
    char *end() {return bytes_.get() + size_;}
    std::size_t capacity() const {return capacity_;}
    std::size_t size() const {return size_;}
    void clear(){size_ = 0;}
    void reserve(std::size_t capacityNew)
    {
        if(capacityNew > capacity_){
            std::unique_ptr<char[]> bytesNew(new char[capacityNew]);
            if(size_ > 0){
                std::memcpy(bytesNew.get(), bytes_.get(), size_);
            }
            bytes_ = std::move(bytesNew);
            capacity_ = capacityNew;
        }
    }
    void resize(std::size_t size)
    {
        reserve(size);
        size_ = size;
        // do no initialize extended elements
    }
};

bool readGzFile(const char *filename, ByteBuffer &buffer)
{
    buffer.clear();

    struct GzFile
    {
        gzFile gzfile_;
        GzFile(const char *filename) : gzfile_(gzopen(filename, "rb")) {}
        ~GzFile() {if(gzfile_){gzclose(gzfile_);}}
        bool isOpen()const {return gzfile_ != NULL;}
        int read(void *buf, unsigned int len){return gzread(gzfile_, buf, len);}
    };
    GzFile gzfile(filename);
    if(!gzfile.isOpen()){
        error(".gz file cannot open(%s)", filename);
        return false;
    }

    const std::size_t SIZE_TO_READ = 1024*1024;
    for(;;){
        if(buffer.capacity() < buffer.size() + SIZE_TO_READ){
            buffer.reserve(buffer.size() + SIZE_TO_READ);
        }
        const int sizeRead = gzfile.read(buffer.end(), buffer.capacity() - buffer.size());
        if(sizeRead == 0){
            break;
        }
        else if(sizeRead < 0){
            error("gzread error(%d)", sizeRead);
            return false;
        }
        buffer.resize(buffer.size() + sizeRead);
    }
    return true;
}

float getColumnAsFloat(
    const std::vector<std::pair<const char *, const char *>> &cols,
    std::size_t index,
    std::optional<float> defaultValue = std::nullopt)
{
    if(index >= cols.size()){
        throw std::out_of_range("index >= cols.size()");
    }

    const char *begin = cols[index].first;
    const char *end = cols[index].second;
    const float value = std::strtof(begin, const_cast<char **>(&end));
    if(value == 0 && end == begin){
        if(defaultValue){
            return *defaultValue;
        }
        else{
            throw std::runtime_error("invalid float string");
        }
    }
    return value;
}

bool extract(FloatOutputFile &outfile,
             const char *filenameCsvGz,
             ByteBuffer &csvbuf,
             unsigned int &numObjects,
             unsigned int &numInvalidLines)
{
    // Read csv file to buffer
    if(!readGzFile(filenameCsvGz, csvbuf)){
        return false;
    }

    // skip first line
    const char *p = csvbuf.get();
    const char *csvend = csvbuf.end();
    for(; p != csvend && *p != '\n' && *p != '\r'; ++p);
    for(; p != csvend && (*p == '\n' || *p == '\r'); ++p);

    // for each line
    std::vector<std::pair<const char *, const char *>> cols;
    while(p != csvend){
        // find first and last position of columns
        cols.clear();
        const char *colbeg = p;
        for(; p != csvend; ++p){
            if(*p == ','){
                cols.emplace_back(colbeg, p);
                colbeg = p + 1;
            }
            else if(*p == '\n' || *p == '\r'){
                break;
            }
        }
        cols.emplace_back(colbeg, p);
        // skip newline characters
        for(; p != csvend && (*p == '\n' || *p == '\r'); ++p);

        // output only necessary columns
        try{
            const float ra = getColumnAsFloat(cols, col_ra);
            const float dec = getColumnAsFloat(cols, col_dec);
            const float g = getColumnAsFloat(cols, col_phot_g_mean_mag);
            const float bp = getColumnAsFloat(cols, col_phot_bp_mean_mag);
            const float rp = getColumnAsFloat(cols, col_phot_rp_mean_mag);
            // https://arxiv.org/pdf/1008.0815.pdf (2)
            const float cxp = bp - rp;
            const float teff = std::max(std::exp(3.999 - 0.654 * cxp + 0.709 * std::pow(cxp, 2) - 0.316 * std::pow(cxp, 3)), 0.0);

            outfile.write(ra);
            outfile.write(dec);
            outfile.write(g);
            outfile.write(bp);
            outfile.write(rp);
            outfile.write(teff);

            ++numObjects;
        }
        catch(const std::exception &e){
            ++numInvalidLines;
        }
    }
    return true;
}


std::string readFirstLine(const char *filename, std::string defaultValue)
{
    char buffer[1024] = "";
    if(std::FILE *file = std::fopen(filename, "r")){
        const std::size_t size = std::fread(buffer, 1, sizeof(buffer), file);
        std::fclose(file);

        for(std::size_t i = 0; i < size; ++i){
            if(buffer[i] == '\n' || buffer[i] == '\r'){
                buffer[i] = '\0';
                return buffer;
            }
        }
        if(size < sizeof(buffer)){
            buffer[size] = '\0';
            return buffer;
        }
        //too long first line of file
    }
    return defaultValue;
}

int main(int argc, char *argv[])
{
    const char * const OUTPUT_FILE = "gaia_ra_dec_g_bp_rp_teff.dat";
    FloatOutputFile outfile(OUTPUT_FILE, 100000);
    if(!outfile.isOpen()){
        error("output file cannot open(%s)", OUTPUT_FILE);
        return false;
    }

    ByteBuffer csvbuf;
    unsigned int numObjects = 0;
    unsigned int numInvalidLines = 0;

    unsigned int numFiles = 0;
    const std::string SOURCE_DIR = readFirstLine("ARCHIVE_DIR", ".") + "/cdn.gea.esac.esa.int/Gaia/gedr3/gaia_source";
    printf("ARCHIVE_DIR=%s\n", SOURCE_DIR.c_str());

    for(auto &p : std::filesystem::directory_iterator(SOURCE_DIR)){
        std::printf("%u\n", numFiles);
        //"GaiaSource_1000172165251650944_1000424567594791808.csv.gz"
        extract(outfile, p.path().c_str(), csvbuf, numObjects, numInvalidLines);
        ++numFiles;
    }

    std::printf("valid:%u invalid:%u error:%u", numObjects, numInvalidLines, numErrors);

}


// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Wei WU (wu.wei.iop@gmail.com)
//
// Template specialization of model parameters for bilayer two-orbital (BLTO) Hubbard model.

template <typename PointGroup>
class ModelParameters<models::TightBindingModel<models::BLTO_lattice<PointGroup>>> {
public:
  ModelParameters() : txp_(0.), tz_(0.),tzp_(0.), t_perp_(0.), Vxz_(0.),delta_(0.) ,U_(0.), V_(0.), V_prime_(0.) {}

  template <typename Concurrency>
  int getBufferSize(const Concurrency& concurrency) const;
  template <typename Concurrency>
  void pack(const Concurrency& concurrency, char* buffer, int buffer_size, int& position) const;
  template <typename Concurrency>
  void unpack(const Concurrency& concurrency, char* buffer, int buffer_size, int& position);

  template <typename ReaderOrWriter>
  void readWrite(ReaderOrWriter& reader_or_writer);

  double get_txp() const {
    return txp_;
  }
  void set_txp(const double t) {
    txp_ = t;
  }

  double get_tz() const {
    return tz_;
  }
  void set_tz(const double t_prime) {
    tz_ = t_prime;
  }

  double get_tzp() const {
    return tzp_;
  }
  void set_tzp(const double t_prime) {
    tzp_ = t_prime;
  }


  double get_t_perp() const {
    return t_perp_;
  }
  void set_t_perp(const double t_perp) {
    t_perp_ = t_perp;
  }
  
  double get_Vxz() const {
    return Vxz_;
  }
  void set_Vxz(const double vv) {
    Vxz_ = vv;
  }
  
   
  double get_delta() const {
    return delta_;
  } 
  
  void set_delta(const double vv) {
    delta_ = vv;
  }  

  double get_U() const {
    return U_;
  }
  void set_U(const double U) {
    U_ = U;
  }

  double get_V() const {
    return V_;
  }
  void set_V(const double V) {
    V_ = V;
  }

  double get_V_prime() const {
    return V_prime_;
  }
  void set_V_prime(const double V_prime) {
    V_prime_ = V_prime;
  }

private:
  double txp_;
  double tz_;
  double tzp_;
  double t_perp_;
  double Vxz_;
  double delta_;
  double U_;
  double V_;
  double V_prime_;

};

template <typename PointGroup>
template <typename Concurrency>
int ModelParameters<models::TightBindingModel<models::BLTO_lattice<PointGroup>>>::getBufferSize(
    const Concurrency& concurrency) const {
  int buffer_size = 0;

  buffer_size += concurrency.get_buffer_size(txp_);
  buffer_size += concurrency.get_buffer_size(tz_);
  buffer_size += concurrency.get_buffer_size(tzp_);
  buffer_size += concurrency.get_buffer_size(Vxz_);
  buffer_size += concurrency.get_buffer_size(t_perp_);
  buffer_size += concurrency.get_buffer_size(delta_);
  buffer_size += concurrency.get_buffer_size(U_);
  buffer_size += concurrency.get_buffer_size(V_);
  buffer_size += concurrency.get_buffer_size(V_prime_);

  return buffer_size;
}

template <typename PointGroup>
template <typename Concurrency>
void ModelParameters<models::TightBindingModel<models::BLTO_lattice<PointGroup>>>::pack(
    const Concurrency& concurrency, char* buffer, int buffer_size, int& position) const {
  concurrency.pack(buffer, buffer_size, position, txp_);
  concurrency.pack(buffer, buffer_size, position, tz_);
  concurrency.pack(buffer, buffer_size, position, tzp_);
  concurrency.pack(buffer, buffer_size, position, Vxz_);
  concurrency.pack(buffer, buffer_size, position, t_perp_);
  concurrency.pack(buffer, buffer_size, position, delta_);
  concurrency.pack(buffer, buffer_size, position, U_);
  concurrency.pack(buffer, buffer_size, position, V_);
  concurrency.pack(buffer, buffer_size, position, V_prime_);
}

template <typename PointGroup>
template <typename Concurrency>
void ModelParameters<models::TightBindingModel<models::BLTO_lattice<PointGroup>>>::unpack(
    const Concurrency& concurrency, char* buffer, int buffer_size, int& position) {
  concurrency.unpack(buffer, buffer_size, position, txp_);
  concurrency.unpack(buffer, buffer_size, position, tz_);
  concurrency.unpack(buffer, buffer_size, position, tzp_);
  concurrency.unpack(buffer, buffer_size, position, Vxz_);
  concurrency.unpack(buffer, buffer_size, position, t_perp_);
  concurrency.unpack(buffer, buffer_size, position, delta_);
  concurrency.unpack(buffer, buffer_size, position, U_);
  concurrency.unpack(buffer, buffer_size, position, V_);
  concurrency.unpack(buffer, buffer_size, position, V_prime_);
}

template <typename PointGroup>
template <typename ReaderOrWriter>
void ModelParameters<models::TightBindingModel<models::BLTO_lattice<PointGroup>>>::readWrite(
    ReaderOrWriter& reader_or_writer) {
  try {
    reader_or_writer.open_group("BLTO-Hubbard-model");

    try {
      reader_or_writer.execute("txp", txp_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("tz", tz_);
    }
    catch (const std::exception& r_e) {
    }
   try {
      reader_or_writer.execute("tzp", tzp_);
    }
    catch (const std::exception& r_e) {
    }
    try{
      reader_or_writer.execute("Vxz", Vxz_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("t-perp", t_perp_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("delta", delta_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("U", U_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("V", V_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("V-prime", V_prime_);
    }
    catch (const std::exception& r_e) {
    }

    reader_or_writer.close_group();
  }
  catch (const std::exception& r_e) {
  }
}

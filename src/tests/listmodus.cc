/*
 *
 *    Copyright (c) 2017-2020,2022-2025
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "vir/test.h"  // This include has to be first

#include "smash/listmodus.h"

#include <cmath>
#include <filesystem>
#include <iostream>
#include <sstream>
#include <string>

#include "setup.h"
#include "smash/oscaroutput.h"
#include "smash/particles.h"

using namespace smash;
static const double accuracy = 5.e-4;
static const std::filesystem::path testoutputpath =
    std::filesystem::absolute(SMASH_TEST_OUTPUT_PATH);
static const auto parameters = Test::default_parameters();

// The type registry is process-global and can only be initialized once.
// Each formation-time test must also work when selected with --only.
static void ensure_particle_types() {
  [[maybe_unused]] static const bool initialized = [] {
    Test::create_stable_smashon_particletypes();
    return true;
  }();
}

static std::filesystem::path create_particlefile(
    const OutputParameters out_par, const int file_number,
    std::vector<ParticleList> &init_particle_vec,
    const int particles_per_event = 10, const int n_events = 1) {
  std::unique_ptr<OutputInterface> osc2013final =
      create_oscar_output("Oscar2013", "Particles", testoutputpath, out_par);
  VERIFY(bool(osc2013final));

  const std::filesystem::path outputfilename = "particle_lists.oscar";
  const std::filesystem::path outputfilepath = testoutputpath / outputfilename;

  // Create random particles
  for (int event = 0; event < n_events; event++) {
    Particles particles;
    for (int i = 0; i < particles_per_event; i++) {
      particles.insert(Test::smashon_random());
    }

    /*
    std::cout << "Initial particles:" << std::endl;
    for (const auto &p : particles) {
      std::cout << p << std::endl;
    }
    */

    init_particle_vec.push_back(particles.copy_to_vector());

    // Print them to file in OSCAR 2013 format
    const double impact_parameter = 2.34;  // just a dummy value here
    const bool empty_event = false;        // just a dummy value as well
    EventInfo default_event_info =
        Test::default_event_info(impact_parameter, empty_event);
    osc2013final->at_eventend(particles, {event, 0}, default_event_info);
  }

  // release and let destructor rename the file
  osc2013final.reset();

  VERIFY(std::filesystem::exists(outputfilepath));
  // Rename the oscar file to match listmodus format
  const std::filesystem::path inputfilepath =
      testoutputpath / ("event" + std::to_string(file_number));
  std::filesystem::rename(outputfilepath, inputfilepath);

  VERIFY(std::filesystem::exists(inputfilepath));

  return inputfilepath;
}

static void create_non_oscar_particlefile(
    const int file_number, std::vector<ParticleList> &init_particle_vec) {
  // Write oscar output, but write only one event and remove all comment lines.
  // This mimics the output of some hydro codes
  OutputParameters out_par = OutputParameters();
  out_par.part_only_final = OutputOnlyFinal::Yes;
  out_par.part_extended = false;
  auto input_path =
      create_particlefile(out_par, file_number, init_particle_vec);

  std::fstream file{input_path};
  const std::filesystem::path tmp_file_str = "tmp_file";
  const std::filesystem::path tmp_path = testoutputpath / tmp_file_str;

  std::ofstream tmp_file;
  tmp_file.open(tmp_path);
  VERIFY(tmp_file.good());
  std::string line;
  while (getline(file, line)) {
    if (line[0] != '#') {
      tmp_file << line << '\n';
    }
  }
  std::filesystem::rename(tmp_path, input_path);
}

static void create_particle_file_with_multiple_particles_at_same_position(
    const OutputParameters out_par) {
  std::unique_ptr<OutputInterface> oscar_output =
      create_oscar_output("Oscar2013", "Particles", testoutputpath, out_par);

  const std::filesystem::path outputfilename = "particle_lists.oscar";
  const std::filesystem::path outputfilepath = testoutputpath / outputfilename;

  // Create some random particles with many at same 4-position
  for (int event = 0; event < 2; event++) {
    Particles particles;
    for (int i = 0; i < 3; i++) {
      particles.insert(Test::smashon_random());
    }
    auto particle = Test::smashon_random();
    for (int i = 0; i < event * 3; i++) {
      particles.insert(particle);
      particle.boost_momentum({0.1, 0.2, 0.3});
    }
    particle = Test::smashon_random();
    for (int i = 0; i < event + 3; i++) {
      particles.insert(particle);
      particle.boost_momentum({0.01, 0.02, 0.03});
    }

    // Print them to file in OSCAR 2013 format
    const double impact_parameter = 2.34;  // just a dummy value here
    const bool empty_event = false;        // just a dummy value as well
    EventInfo default_event_info =
        Test::default_event_info(impact_parameter, empty_event);
    oscar_output->at_eventend(particles, {event, 0}, default_event_info);
  }

  // release and let destructor rename the file
  oscar_output.reset();

  VERIFY(std::filesystem::exists(outputfilepath));
  // Rename the oscar file to match listmodus format
  const std::filesystem::path inputfilepath = testoutputpath / "event0";
  std::filesystem::rename(outputfilepath, inputfilepath);
}

static ListModus create_list_modus_for_test() {
  Configuration config{R"(
    Modi:
      List:
        File_Directory: ToBeSet
        File_Prefix: event
    )"};
  config.set_value(InputKeys::modi_list_fileDirectory, testoutputpath.string());
  return ListModus(std::move(config), parameters);
}

static ListBoxModus create_list_box_modus_for_test() {
  Configuration config{R"(
    Modi:
      ListBox:
        File_Directory: ToBeSet
        File_Prefix: event
        Length: 3
    )"};
  config.set_value(InputKeys::modi_listBox_fileDirectory,
                   testoutputpath.string());
  return ListBoxModus(std::move(config), parameters);
}

static ListModus create_list_modus_with_single_file_for_test() {
  Configuration config{R"(
    Modi:
      List:
        File_Directory: ToBeSet
        Filename: event0
    )"};
  config.set_value(InputKeys::modi_list_fileDirectory, testoutputpath.string());
  return ListModus(std::move(config), parameters);
}

TEST(directory_is_created) {
  std::filesystem::create_directories(testoutputpath);
  VERIFY(std::filesystem::exists(testoutputpath));
}

TEST(create_particle_types) { ensure_particle_types(); }

static void compare_fourvector(const FourVector &a, const FourVector &b) {
  COMPARE_ABSOLUTE_ERROR(a.x0(), b.x0(), accuracy);
  COMPARE_ABSOLUTE_ERROR(a.x1(), b.x1(), accuracy);
  COMPARE_ABSOLUTE_ERROR(a.x2(), b.x2(), accuracy);
  COMPARE_ABSOLUTE_ERROR(a.x3(), b.x3(), accuracy);
}

TEST(list_from_non_oscar_output) {
  std::vector<ParticleList> init_particles;
  create_non_oscar_particlefile(0, init_particles);
  ListModus list_modus = create_list_modus_for_test();

  // Read the file with list modus
  Particles particles_read;
  list_modus.initial_conditions(&particles_read, parameters);

  /*
  std::cout << "Particles from list modus:" << std::endl;
  for (const auto &p : particles_read) {
    std::cout << p << std::endl;
  }
  */

  // Scroll particles back to the earliest time, as list modus should do
  double earliest_t = 1.e8;
  for (const auto &particle : init_particles[0]) {
    if (particle.position().x0() < earliest_t) {
      earliest_t = particle.position().x0();
    }
  }
  for (auto &particle : init_particles[0]) {
    const double t = particle.position().x0();
    const FourVector u(1.0, particle.velocity());
    particle.set_formation_time(t);
    particle.set_4position(particle.position() + u * (earliest_t - t));
  }

  COMPARE(particles_read.size(), init_particles[0].size());
  ParticleList p_init = init_particles[0];
  ParticleList p_fin = particles_read.copy_to_vector();
  for (size_t i = 0; i < p_fin.size(); i++) {
    ParticleData a = p_init.back();
    ParticleData b = p_fin.back();
    p_init.pop_back();
    p_fin.pop_back();
    compare_fourvector(a.momentum(), b.momentum());
    compare_fourvector(a.position(), b.position());
    COMPARE(a.id(), b.id());
    COMPARE_ABSOLUTE_ERROR(a.formation_time(), b.formation_time(), accuracy);
    COMPARE(a.pdgcode(), b.pdgcode());
  }
}

TEST(listbox_creation_from_non_oscar_output) {
  std::vector<ParticleList> init_particles;
  create_non_oscar_particlefile(0, init_particles);
  ListBoxModus list_box_modus = create_list_box_modus_for_test();

  // Read the file with list modus
  Particles particles_read;
  list_box_modus.initial_conditions(&particles_read, parameters);
}

TEST(multiple_file_non_oscar_output) {
  ListModus list_modus = create_list_modus_for_test();

  std::vector<ParticleList> init_particles;

  constexpr size_t max_events = 10;
  for (size_t i = 0; i < max_events; i++) {
    create_non_oscar_particlefile(i, init_particles);
  }
  COMPARE(init_particles.size(), max_events);

  // Read particles with list modus
  for (size_t current_event = 0; current_event < max_events; current_event++) {
    Particles particles_read;
    list_modus.initial_conditions(&particles_read, parameters);

    /*
    std::cout << "Particles from list modus:" << std::endl;
    for (const auto &p : particles_read) {
      std::cout << p << std::endl;
    }
    */

    // Scroll particles back to the earliest time, as list modus should do
    double earliest_t = 1.e8;
    for (const auto &particle : init_particles[current_event]) {
      if (particle.position().x0() < earliest_t) {
        earliest_t = particle.position().x0();
      }
    }
    for (auto &particle : init_particles[current_event]) {
      const double t = particle.position().x0();
      const FourVector u(1.0, particle.velocity());
      particle.set_formation_time(t);
      particle.set_4position(particle.position() + u * (earliest_t - t));
    }

    COMPARE(particles_read.size(), init_particles[current_event].size());
    ParticleList p_init = init_particles[current_event];
    ParticleList p_fin = particles_read.copy_to_vector();
    for (size_t i = 0; i < p_fin.size(); i++) {
      ParticleData a = p_init.back();
      ParticleData b = p_fin.back();
      p_init.pop_back();
      p_fin.pop_back();
      compare_fourvector(a.momentum(), b.momentum());
      compare_fourvector(a.position(), b.position());
      COMPARE(a.id(), b.id());
      COMPARE_ABSOLUTE_ERROR(a.formation_time(), b.formation_time(), accuracy);
      COMPARE(a.pdgcode(), b.pdgcode());
    }
  }
}

TEST(list_from_oscar2013_output) {
  // Create OSCAR 2013 output
  OutputParameters out_par = OutputParameters();
  out_par.part_only_final = OutputOnlyFinal::Yes;
  out_par.part_extended = false;
  std::vector<ParticleList> init_particles;
  create_particlefile(out_par, 0, init_particles, 10, 1);
  ListModus list_modus = create_list_modus_for_test();

  // Read the file with list modus
  Particles particles_read;
  list_modus.initial_conditions(&particles_read, parameters);

  /*
  std::cout << "Particles from list modus:" << std::endl;
  for (const auto &p : particles_read) {
    std::cout << p << std::endl;
  }
  */

  // Scroll particles back to the earliest time, as list modus should do
  double earliest_t = 1.e8;
  for (const auto &particle : init_particles[0]) {
    if (particle.position().x0() < earliest_t) {
      earliest_t = particle.position().x0();
    }
  }
  for (auto &particle : init_particles[0]) {
    const double t = particle.position().x0();
    const FourVector u(1.0, particle.velocity());
    particle.set_formation_time(t);
    particle.set_4position(particle.position() + u * (earliest_t - t));
  }

  COMPARE(particles_read.size(), init_particles[0].size());
  ParticleList p_init = init_particles[0];
  ParticleList p_fin = particles_read.copy_to_vector();
  for (size_t i = 0; i < p_fin.size(); i++) {
    ParticleData a = p_init.back();
    ParticleData b = p_fin.back();
    p_init.pop_back();
    p_fin.pop_back();
    compare_fourvector(a.momentum(), b.momentum());
    compare_fourvector(a.position(), b.position());
    COMPARE(a.id(), b.id());
    COMPARE_ABSOLUTE_ERROR(a.formation_time(), b.formation_time(), accuracy);
    COMPARE(a.pdgcode(), b.pdgcode());
  }
}

TEST(multiple_files_one_event) {
  OutputParameters out_par = OutputParameters();
  out_par.part_only_final = OutputOnlyFinal::Yes;
  out_par.part_extended = false;
  std::vector<ParticleList> init_particles;
  constexpr int events_per_file = 1;
  constexpr int particles_per_event = 2;
  constexpr int n_files = 2;
  for (int i = 0; i < n_files; i++) {
    create_particlefile(out_par, i, init_particles, particles_per_event,
                        events_per_file);
  }
  ListModus list_modus = create_list_modus_for_test();

  for (int i = 0; i < events_per_file * n_files; i++) {
    Particles particles_read;
    list_modus.initial_conditions(&particles_read, parameters);

    // Scroll particles back to the earliest time, as list modus should do
    double earliest_t = 1.e8;
    for (const auto &particle : init_particles[i]) {
      if (particle.position().x0() < earliest_t) {
        earliest_t = particle.position().x0();
      }
    }
    for (auto &particle : init_particles[i]) {
      const double t = particle.position().x0();
      const FourVector u(1.0, particle.velocity());
      particle.set_formation_time(t);
      particle.set_4position(particle.position() + u * (earliest_t - t));
    }

    COMPARE(particles_read.size(), init_particles[i].size());
    ParticleList p_init = init_particles[i];
    ParticleList p_fin = particles_read.copy_to_vector();
    for (size_t j = 0; j < p_fin.size(); j++) {
      ParticleData a = p_init.back();
      ParticleData b = p_fin.back();
      p_init.pop_back();
      p_fin.pop_back();
      compare_fourvector(a.momentum(), b.momentum());
      compare_fourvector(a.position(), b.position());
      COMPARE(a.id(), b.id());
      COMPARE_ABSOLUTE_ERROR(a.formation_time(), b.formation_time(), accuracy);
      COMPARE(a.pdgcode(), b.pdgcode());
    }
  }
}

TEST(multiple_files_multiple_events) {
  OutputParameters out_par = OutputParameters();
  out_par.part_only_final = OutputOnlyFinal::Yes;
  out_par.part_extended = false;

  std::vector<ParticleList> init_particles;
  constexpr int events_per_file = 5;
  constexpr int particles_per_event = 10;
  constexpr int n_files = 5;

  for (int i = 0; i < n_files; i++) {
    create_particlefile(out_par, i, init_particles, particles_per_event,
                        events_per_file);
  }
  ListModus list_modus = create_list_modus_for_test();

  for (int i = 0; i < events_per_file * n_files; i++) {
    Particles particles_read;
    list_modus.initial_conditions(&particles_read, parameters);

    // Scroll particles back to the earliest time, as ListModus should do
    double earliest_t = 1.e8;
    for (const auto &particle : init_particles[i]) {
      if (particle.position().x0() < earliest_t) {
        earliest_t = particle.position().x0();
      }
    }
    for (auto &particle : init_particles[i]) {
      const double t = particle.position().x0();
      const FourVector u(1.0, particle.velocity());
      particle.set_formation_time(t);
      particle.set_4position(particle.position() + u * (earliest_t - t));
    }

    COMPARE(particles_read.size(), init_particles[i].size());
    ParticleList p_init = init_particles[i];
    ParticleList p_fin = particles_read.copy_to_vector();
    for (size_t j = 0; j < p_fin.size(); j++) {
      ParticleData a = p_init.back();
      ParticleData b = p_fin.back();
      p_init.pop_back();
      p_fin.pop_back();
      compare_fourvector(a.momentum(), b.momentum());
      compare_fourvector(a.position(), b.position());
      COMPARE(a.id(), b.id());
      COMPARE_ABSOLUTE_ERROR(a.formation_time(), b.formation_time(), accuracy);
      COMPARE(a.pdgcode(), b.pdgcode());
    }
  }
}

TEST(multiple_events_in_file) {
  OutputParameters out_par = OutputParameters();
  out_par.part_only_final = OutputOnlyFinal::Yes;
  out_par.part_extended = false;
  constexpr int max_events = 2;
  constexpr int particles_per_event = 10;
  std::vector<ParticleList> init_particles;
  create_particlefile(out_par, 0, init_particles, particles_per_event,
                      max_events);
  ListModus list_modus = create_list_modus_with_single_file_for_test();

  for (int cur_event = 0; cur_event < max_events; cur_event++) {
    // Read the file with list modus
    Particles particles_read;
    list_modus.initial_conditions(&particles_read, parameters);

    /*
    std::cout << "Particles from list modus:" << std::endl;
    for (const auto &p : particles_read) {
      std::cout << p << std::endl;
    }
    */

    // Scroll particles back to the earliest time, as list modus should do
    double earliest_t = 1.e8;
    for (const auto &particle : init_particles[cur_event]) {
      if (particle.position().x0() < earliest_t) {
        earliest_t = particle.position().x0();
      }
    }
    for (auto &particle : init_particles[cur_event]) {
      const double t = particle.position().x0();
      const FourVector u(1.0, particle.velocity());
      particle.set_formation_time(t);
      particle.set_4position(particle.position() + u * (earliest_t - t));
    }

    COMPARE(particles_read.size(), init_particles[cur_event].size());
    ParticleList p_init = init_particles[cur_event];
    ParticleList p_fin = particles_read.copy_to_vector();
    for (size_t i = 0; i < p_fin.size(); i++) {
      ParticleData a = p_init.back();
      ParticleData b = p_fin.back();
      p_init.pop_back();
      p_fin.pop_back();
      compare_fourvector(a.momentum(), b.momentum());
      compare_fourvector(a.position(), b.position());
      COMPARE(a.id(), b.id());
      COMPARE_ABSOLUTE_ERROR(a.formation_time(), b.formation_time(), accuracy);
      COMPARE(a.pdgcode(), b.pdgcode());
    }
  }
}

TEST(try_create_particle_func) {
  ListModus list_modus = create_list_modus_for_test();
  Particles particles;
  ParticleList plist_init, plist_fin;
  const int npart = 10;
  const double m0 = Test::smashon_mass;

  for (int i = 0; i < npart; i++) {
    ParticleData smashon = Test::smashon_random();
    plist_init.push_back(smashon);
    FourVector r = smashon.position(), p = smashon.momentum();
    PdgCode pdg = smashon.pdgcode();
    list_modus.try_create_particle(particles, pdg, r.x0(), r.x1(), r.x2(),
                                   r.x3(), m0, p.x0(), p.x1(), p.x2(), p.x3());
  }
  plist_fin = particles.copy_to_vector();
  for (int i = 0; i < npart; i++) {
    ParticleData a = plist_init.back();
    ParticleData b = plist_fin.back();
    plist_init.pop_back();
    plist_fin.pop_back();
    compare_fourvector(a.momentum(), b.momentum());
    compare_fourvector(a.position(), b.position());
    COMPARE_ABSOLUTE_ERROR(b.formation_time(), a.position().x0(), accuracy);
    COMPARE(a.pdgcode(), b.pdgcode());
  }

  // Create stable smashons some with a mass discrepancy and some off-shell
  // because of their energy. Test if they are returned on-shell with pole mass.
  particles.reset();
  plist_init.clear();
  plist_fin.clear();
  for (int i = 0; i < 2 * npart; i++) {
    ParticleData smashon = Test::smashon_random();
    plist_init.push_back(smashon);
    FourVector r = smashon.position(), p = smashon.momentum();
    PdgCode pdg = smashon.pdgcode();
    const auto creation_mass = m0 + (i < npart);
    const auto creation_energy = p.x0() + (i >= npart);
    list_modus.try_create_particle(particles, pdg, r.x0(), r.x1(), r.x2(),
                                   r.x3(), creation_mass, creation_energy,
                                   p.x1(), p.x2(), p.x3());
  }
  plist_fin = particles.copy_to_vector();
  for (int i = 0; i < npart; i++) {
    ParticleData a = plist_init.back();
    ParticleData b = plist_fin.back();
    plist_init.pop_back();
    plist_fin.pop_back();
    // Test smashon should be on mass shell with pole mass
    COMPARE_ABSOLUTE_ERROR(a.momentum().abs(), m0, accuracy);
    // The smashon read from the list should have the pole mass,
    // but still obey E^2 - p^2 = m^2.
    COMPARE_ABSOLUTE_ERROR(b.momentum().abs(), m0, accuracy);
    compare_fourvector(a.position(), b.position());
    COMPARE_ABSOLUTE_ERROR(b.formation_time(), a.position().x0(), accuracy);
    COMPARE(a.pdgcode(), b.pdgcode());
  }
}

TEST(try_create_particle_with_spin_func) {
  // Create a ListModus with spin interactions enabled
  const ExperimentParameters parameters_with_spin{
      std::make_unique<UniformClock>(0., 0.1, 300.0),
      std::make_unique<UniformClock>(0., 1., 300.0),
      1,  // ensembles
      1,  // testparticles
      DerivativesMode::CovariantGaussian,
      RestFrameDensityDerivativesMode::Off,
      FieldDerivativesMode::ChainRule,
      SmearingMode::CovariantGaussian,
      1.0,
      4.0,
      0.333333,
      2.0,
      CollisionCriterion::Geometric,
      true,
      Test::all_reactions_included(),
      Test::no_multiparticle_reactions(),
      false,  // strings
      1.0,
      NNbarTreatment::NoAnnihilation,
      0.,
      false,
      -1.0,
      200.0,
      2.5,
      1.0,
      false,
      false,
      true,
      SpinInteractionType::On,
      std::nullopt};

  Configuration config{R"(
    Modi:
      List:
        File_Directory: ToBeSet
        File_Prefix: event
    )"};
  config.set_value(InputKeys::modi_list_fileDirectory, testoutputpath.string());
  config.set_value(
      InputKeys::modi_list_optionalQuantities,
      std::vector<std::string>{"spin0", "spinx", "spiny", "spinz"});
  ListModus list_modus = ListModus(std::move(config), parameters_with_spin);
  Particles particles;
  const ParticleData smashon = Test::smashon_random();
  const FourVector r = smashon.position(), p = smashon.momentum();
  const PdgCode pdg = smashon.pdgcode();
  // Spin vector components as 4-vector and optional values
  const FourVector spin_vec(0.1, 0.2, 0.3, 0.4);
  const std::vector<std::string> opt_vals = {"0.1", "0.2", "0.3", "0.4"};
  list_modus.try_create_particle(particles, pdg, r.x0(), r.x1(), r.x2(), r.x3(),
                                 Test::smashon_mass, p.x0(), p.x1(), p.x2(),
                                 p.x3(), opt_vals);

  const ParticleList plist = particles.copy_to_vector();
  const ParticleData created = plist.back();
  compare_fourvector(created.spin_vector(), spin_vec);
}

TEST(imported_production_and_formation_times) {
  ensure_particle_types();
  const auto input_path = testoutputpath / "formation_times";
  std::filesystem::create_directories(input_path);
  for (const bool reverse_fields : {false, true}) {
    Configuration config{R"(
      Modi:
        List:
          File_Directory: ToBeSet
          File_Prefix: event
      )"};
    config.set_value(InputKeys::modi_list_fileDirectory, input_path.string());
    config.set_value(InputKeys::modi_list_optionalQuantities,
                     reverse_fields
                         ? std::vector<std::string>{"form_time", "prod_time"}
                         : std::vector<std::string>{"prod_time", "form_time"});
    if (reverse_fields) {
      config.set_value(
          InputKeys::collTerm_xsecFormationScaling_formationTime_mode,
          std::string{"OSCAR"});
      config.set_value(
          InputKeys::collTerm_xsecFormationScaling_formationTime_properTime,
          100.0);
    }
    ListModus list_modus(std::move(config), parameters);
    Particles particles;
    const double mass = Test::smashon_mass;
    const double px = 0.4;
    const double energy = std::sqrt(mass * mass + px * px);
    for (int i = 0; i < 3; ++i) {
      const double prod_time = i == 2 ? 0.0 : 2.0 + i;
      const double form_time = i == 2 ? 0.0 : 9.0 + 2.0 * i;
      const std::vector<std::string> times =
          reverse_fields
              ? std::vector<std::string>{std::to_string(form_time),
                                         std::to_string(prod_time)}
              : std::vector<std::string>{std::to_string(prod_time),
                                         std::to_string(form_time)};
      list_modus.try_create_particle(particles, PdgCode("661"), 0.0, i, 0.0, 0.0,
                                     mass, energy, px, 0.0, 0.0, times);
      const ParticleData created = particles.copy_to_vector().back();
      COMPARE(created.begin_formation_time(), prod_time);
      COMPARE(created.formation_time(), form_time);
      COMPARE(created.position().x0(), 0.0);
    }
  }
}

TEST(constant_proper_formation_times) {
  ensure_particle_types();
  const auto input_path = testoutputpath / "constant_proper_times";
  std::filesystem::create_directories(input_path);
  for (const bool include_form_time : {false, true}) {
    for (const double proper_time : {0.0, 1.6}) {
      Configuration config{R"(
        Modi:
          List:
            File_Directory: ToBeSet
            File_Prefix: event
        Collision_Term:
          Cross_Section_Formation_Scaling:
            Formation_Time:
              Mode: Constant_Proper_Time
        )"};
      config.set_value(InputKeys::modi_list_fileDirectory, input_path.string());
      config.set_value(
          InputKeys::collTerm_xsecFormationScaling_formationTime_properTime,
          proper_time);
      std::vector<std::string> fields{"prod_time_status", "xsecfac", "prod_time"};
      if (include_form_time) fields.push_back("form_time");
      config.set_value(InputKeys::modi_list_optionalQuantities, fields);
      ListModus list_modus(std::move(config), parameters);
      Particles particles;
      for (const double gamma : {1.0, 1.25, 5.0}) {
        const double mass = Test::smashon_mass;
        const double px = mass * std::sqrt(gamma * gamma - 1.0);
        std::vector<std::string> values{"1", "0", "2.5"};
        if (include_form_time) values.push_back("999");
        list_modus.try_create_particle(
            particles, PdgCode("661"), 0.0, gamma, 0.0, 0.0, mass,
            gamma * mass, px, 0.0, 0.0, values);
        const auto created = particles.copy_to_vector().back();
        COMPARE(created.begin_formation_time(), 2.5);
        COMPARE_ABSOLUTE_ERROR(created.formation_time(),
                               2.5 + gamma * proper_time, 1.e-12);
        COMPARE(created.position().x0(), 0.0);
        COMPARE(created.position().x1(), gamma);
        COMPARE(created.initial_xsec_scaling_factor(), 0.0);
        COMPARE(created.xsec_scaling_factor(created.formation_time()), 1.0);
      }
      // Spectators must not acquire a new formation interval, even if the
      // imported endpoint says otherwise. Also verify available zero times.
      for (const int status : {0, 1, 2}) {
        std::vector<std::string> values{
            std::to_string(status), "1", status == 0 ? "7" : "0"};
        if (include_form_time) values.push_back("999");
        list_modus.try_create_particle(
            particles, PdgCode("661"), 0.0, status, 0.0, 0.0,
            Test::smashon_mass, Test::smashon_mass, 0.0, 0.0, 0.0, values);
        const auto created = particles.copy_to_vector().back();
        COMPARE(created.begin_formation_time(), 0.0);
        COMPARE(created.formation_time(), status == 2 ? 0.0 : proper_time);
        COMPARE(created.xsec_scaling_factor(), 1.0);
      }
    }
  }
}

TEST(constant_proper_time_listbox_file) {
  ensure_particle_types();
  const auto input_path = testoutputpath / "constant_proper_listbox";
  std::filesystem::create_directories(input_path);
  {
    std::ofstream file(input_path / "event0");
    file << "# event 0 start\n"
         << "0 1 0 0 0.123 0.15375 0.09225 0 0 661 2.5 0 1\n"
         << "0 2 0 0 0.123 0.123 0 0 0 661 0 1 2\n"
         << "# event 0 end\n";
  }
  Configuration config{R"(
    Modi:
      ListBox:
        File_Directory: ToBeSet
        Filename: event0
        Length: 10.0
        Optional_Quantities: [prod_time, xsecfac, prod_time_status]
    Collision_Term:
      Cross_Section_Formation_Scaling:
        Formation_Time:
          Mode: Constant_Proper_Time
          Proper_Time: 1.6
    )"};
  config.set_value(InputKeys::modi_listBox_fileDirectory, input_path.string());
  ListBoxModus list_modus(std::move(config), parameters);
  // Exercise file parsing, startup validation, and single-file replay.
  for (int replay = 0; replay < 2; ++replay) {
    Particles particles;
    list_modus.initial_conditions(&particles, parameters);
    const auto created = particles.copy_to_vector();
    COMPARE(created.size(), 2u);
    COMPARE_ABSOLUTE_ERROR(created[0].formation_time(), 4.5, 1.e-12);
    COMPARE(created[1].formation_time(), 0.0);
  }
}

TEST(invalid_formation_time_configuration) {
  const auto input_path = testoutputpath / "invalid_formation_mode";
  std::filesystem::create_directories(input_path);
  for (const std::string &setting : {"{Mode: Unknown}",
                                     "{Mode: Constant_Proper_Time, Proper_Time: -1}",
                                     "{Mode: Constant_Proper_Time, Proper_Time: .nan}",
                                     "{Mode: Constant_Proper_Time, Proper_Time: .inf}"}) {
    const std::string yaml = std::string{
        "Modi: {List: {File_Directory: ToBeSet, File_Prefix: event, "
        "Optional_Quantities: [prod_time, form_time]}}\n"
        "Collision_Term: {Cross_Section_Formation_Scaling: {Formation_Time: "} +
        setting + "}}";
    Configuration config{yaml.c_str()};
    config.set_value(InputKeys::modi_list_fileDirectory, input_path.string());
    bool rejected = false;
    try {
      ListModus list_modus(std::move(config), parameters);
    } catch (const std::invalid_argument &) {
      rejected = true;
    }
    VERIFY(rejected);
  }
}

TEST(constant_proper_time_requires_production) {
  Configuration config{R"(
    Modi:
      List:
        File_Directory: ToBeSet
        File_Prefix: event
        Optional_Quantities: [form_time]
    Collision_Term:
      Cross_Section_Formation_Scaling:
        Formation_Time:
          Mode: Constant_Proper_Time
    )"};
  config.set_value(InputKeys::modi_list_fileDirectory, testoutputpath.string());
  bool rejected = false;
  try {
    ListModus list_modus(std::move(config), parameters);
  } catch (const std::invalid_argument &error) {
    rejected = true;
    VERIFY(std::string(error.what()).find("requires prod_time") !=
           std::string::npos);
  }
  VERIFY(rejected);
}

TEST(imported_production_time_status) {
  ensure_particle_types();
  const auto input_path = testoutputpath / "production_time_status";
  std::filesystem::create_directories(input_path);
  // Restore stderr even if a test assertion or the reader throws.
  struct CaptureStderr {
    std::ostringstream output;
    std::streambuf *saved = std::cerr.rdbuf(output.rdbuf());
    ~CaptureStderr() { std::cerr.rdbuf(saved); }
  };
  for (const bool status_first : {false, true}) {
    {
      std::ofstream file(input_path / "event0");
      file << "#!OSCAR2013 particle_lists t x y z mass p0 px py pz pdg ID charge ";
      file << (status_first ? "prod_time_status prod_time form_time xsecfac\n"
                            : "prod_time form_time xsecfac prod_time_status\n");
      file << "# event 0 out 5\n";
      for (int i = 0; i < 5; ++i) {
        const int status = i < 2 ? 0 : i < 4 ? 1 : 2;
        const double production = i == 1 ? 7.0 : i == 3 ? 3.0 : 0.0;
        const double formation = i == 4 ? 0.0 : 9.0;
        const double factor = i == 4 ? 1.0 : 0.0;
        file << "0 " << i << " 0 0 " << Test::smashon_mass << ' '
             << Test::smashon_mass << " 0 0 0 661 " << i << " 0 ";
        if (status_first) file << status << ' ';
        file << production << ' ' << formation << ' ' << factor;
        if (!status_first) file << ' ' << status;
        file << '\n';
      }
      file << "# event 0 end\n";
    }
    Configuration config{R"(
      Modi:
        List:
          File_Directory: ToBeSet
          File_Prefix: event
      )"};
    config.set_value(InputKeys::modi_list_fileDirectory, input_path.string());
    config.set_value(InputKeys::modi_list_optionalQuantities,
                     status_first
                         ? std::vector<std::string>{"ID", "charge", "prod_time_status",
                                                    "prod_time", "form_time", "xsecfac"}
                         : std::vector<std::string>{"ID", "charge", "prod_time",
                                                    "form_time", "xsecfac", "prod_time_status"});
    CaptureStderr capture;
    ListModus list_modus(std::move(config), parameters);
    Particles particles;
    list_modus.initial_conditions(&particles, parameters);
    const auto created = particles.copy_to_vector();
    COMPARE(created.size(), 5u);
    for (size_t i = 0; i < created.size(); ++i) {
      COMPARE(created[i].begin_formation_time(), i == 3 ? 3.0 : 0.0);
      COMPARE(created[i].formation_time(), i == 4 ? 0.0 : 9.0);
    }
    const std::string warning = "prod_time_status=0); using prod_time=0.";
    const auto text = capture.output.str();
    size_t count = 0;
    for (size_t pos = 0; (pos = text.find(warning, pos)) != std::string::npos;
         pos += warning.size()) {
      ++count;
    }
    // Two missing times warn during startup validation and again on import.
    // Available zero and not-applicable rows produce no warning.
    COMPARE(count, 4u);
    VERIFY(text.find("particle PDG 661") != std::string::npos);
  }
}

TEST(invalid_production_time_status) {
  ensure_particle_types();
  const auto input_path = testoutputpath / "invalid_production_time_status";
  std::filesystem::create_directories(input_path);
  Configuration config{R"(
    Modi:
      List:
        File_Directory: ToBeSet
        File_Prefix: event
    )"};
  config.set_value(InputKeys::modi_list_fileDirectory, input_path.string());
  config.set_value(InputKeys::modi_list_optionalQuantities,
                   std::vector<std::string>{"prod_time", "form_time", "prod_time_status"});
  ListModus list_modus(std::move(config), parameters);
  for (const auto &values : {std::vector<std::string>{"2", "9", "-1"},
                             std::vector<std::string>{"2", "9", "3"},
                             std::vector<std::string>{"2", "9", "1.5"},
                             std::vector<std::string>{"2", "9", "1junk"},
                             std::vector<std::string>{"2", "9", ""},
                             std::vector<std::string>{"2", "9"},
                             std::vector<std::string>{"nan", "9", "1"}}) {
    Particles particles;
    bool rejected = false;
    try {
      list_modus.try_create_particle(particles, PdgCode("661"), 0.0, 0.0, 0.0,
                                     0.0, Test::smashon_mass, Test::smashon_mass,
                                     0.0, 0.0, 0.0, values);
    } catch (const std::invalid_argument &) {
      rejected = true;
    } catch (const std::out_of_range &) {
      rejected = true;
    }
    VERIFY(rejected);
    VERIFY(particles.is_empty());
  }
}

TEST(missing_production_or_formation_fields) {
  for (const auto &fields : {std::vector<std::string>{},
                             std::vector<std::string>{"prod_time"},
                             std::vector<std::string>{"form_time"}}) {
    Configuration config{R"(
      Modi:
        List:
          File_Directory: ToBeSet
          File_Prefix: event
      )"};
    config.set_value(InputKeys::modi_list_fileDirectory, testoutputpath.string());
    config.set_value(InputKeys::modi_list_optionalQuantities, fields);
    bool rejected = false;
    try {
      ListModus list_modus(std::move(config), parameters);
    } catch (const std::invalid_argument &error) {
      rejected = true;
      VERIFY(std::string(error.what()).find("prod_time and form_time") !=
             std::string::npos);
    }
    VERIFY(rejected);
  }
}

TEST(missing_production_or_formation_values) {
  ensure_particle_types();
  const auto input_path = testoutputpath / "formation_times";
  std::filesystem::create_directories(input_path);
  Configuration config{R"(
    Modi:
      List:
        File_Directory: ToBeSet
        File_Prefix: event
    )"};
  config.set_value(InputKeys::modi_list_fileDirectory, input_path.string());
  config.set_value(InputKeys::modi_list_optionalQuantities,
                   std::vector<std::string>{"prod_time", "form_time"});
  ListModus list_modus(std::move(config), parameters);
  for (const auto &values : {std::vector<std::string>{},
                             std::vector<std::string>{"2"},
                             std::vector<std::string>{"", "9"},
                             std::vector<std::string>{"2", ""}}) {
    Particles particles;
    bool rejected = false;
    try {
      list_modus.try_create_particle(particles, PdgCode("661"), 0.0, 0.0, 0.0,
                                     0.0, Test::smashon_mass,
                                     Test::smashon_mass, 0.0, 0.0, 0.0, values);
    } catch (const std::invalid_argument &) {
      rejected = true;
    } catch (const std::out_of_range &) {
      rejected = true;
    }
    VERIFY(rejected);
    VERIFY(particles.is_empty());
  }
}

TEST_CATCH(missing_time_in_particle_file, ListModus::LoadFailure) {
  ensure_particle_types();
  const auto input_path = testoutputpath / "missing_time";
  std::filesystem::create_directories(input_path);
  {
    std::ofstream file(input_path / "event0");
    // The second particle omits form_time; validation must fail at startup.
    file << "# event 0 start\n"
         << "0 0 0 0 0.123 0.123 0 0 0 661 2 9\n"
         << "0 1 0 0 0.123 0.123 0 0 0 661 3\n"
         << "# event 0 end\n";
  }
  Configuration config{R"(
    Modi:
      List:
        File_Directory: ToBeSet
        File_Prefix: event
    )"};
  config.set_value(InputKeys::modi_list_fileDirectory, input_path.string());
  config.set_value(InputKeys::modi_list_optionalQuantities,
                   std::vector<std::string>{"prod_time", "form_time"});
  ListModus list_modus(std::move(config), parameters);
}

TEST_CATCH(insert_optional_quantities_to_func, std::out_of_range) {
  ListModus list_modus = create_list_modus_for_test();
  Particles particles;
  const double m0 = Test::smashon_mass;
  const ParticleData smashon = Test::smashon_random();
  const FourVector r = smashon.position(), p = smashon.momentum();
  const PdgCode pdg = smashon.pdgcode();

  // Create a particle with any optional_quantities handed to the function as
  // last argument other than the default {} to trigger an out_of_range error.
  list_modus.try_create_particle(particles, pdg, r.x0(), r.x1(), r.x2(), r.x3(),
                                 m0, p.x0(), p.x1(), p.x2(), p.x3(), {"ID"});
}

TEST_CATCH(create_particle_with_nan, std::invalid_argument) {
  ListModus list_modus = create_list_modus_for_test();
  Particles particles;
  const double m0 = Test::smashon_mass;
  const ParticleData smashon = Test::smashon_random();
  const FourVector r = smashon.position(), p = smashon.momentum();
  const PdgCode pdg = smashon.pdgcode();

  // Create a particle with either a NAN value in the position
  // to trigger an invalid_argument error.
  list_modus.try_create_particle(particles, pdg, NAN, r.x1(), r.x2(), r.x3(),
                                 m0, p.x0(), p.x1(), p.x2(), p.x3());
}

TEST_CATCH(create_particles_at_same_position, ListModus::InvalidEvents) {
  const OutputParameters out_par = OutputParameters();
  create_particle_file_with_multiple_particles_at_same_position(out_par);
  ListModus list_modus = create_list_modus_with_single_file_for_test();
}

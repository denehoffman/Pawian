#ifndef PARTICLETABLE_HH
#define PARTICLETABLE_HH

#include <map>
#include <string>

class Particle;

typedef std::pair<const std::string, Particle*> pdt_value_type;


class ParticleTable
{
public:

  ParticleTable();
  ~ParticleTable();

  Particle* particle(std::string name);

  bool addParticle(Particle* newParticle);
  bool clone(std::string newName, std::string oldName);

  bool modifyMass(std::string name, double newMass);
  bool modifyWidth(std::string name, double newWidth);

  void print(std::ostream& o) const;

private:
  std::map<const std::string, Particle*> particles;

};

#endif


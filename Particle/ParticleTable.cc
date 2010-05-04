#include "Particle/ParticleTable.hh"
#include "Particle/Particle.hh"

ParticleTable::ParticleTable()
{
}

ParticleTable::~ParticleTable()
{
}

bool ParticleTable::addParticle(Particle* newParticle)
{
  std::string* key = new std::string(newParticle->name());
  particles.insert(pdt_value_type(*key, newParticle));
  return true;
}

Particle* ParticleTable::particle(std::string name)
{
  std::map<const std::string, Particle*>::const_iterator iter = particles.find(name);
  if (iter != particles.end())
    return iter->second;
  else
    return 0;
}


bool ParticleTable::clone(std::string newName, std::string oldName)
{
  std::map<const std::string, Particle*>::const_iterator iter = particles.find(newName);
  if (iter != particles.end()) {
    std::cerr << "ParticleTable: can not clone to already existing particle " << newName << std::endl;
    return false; // cannot clone to an existing particle
  }
  
  Particle* oldParticle = particle(oldName);
  if (0 == oldParticle) {
    std::cerr << "ParticleTable: cannot clone, " << oldName << " does not exist!" << std::endl;
    return false;
  }

  ParticleData* newData = new ParticleData(oldParticle->data());
  newData->name = newName;
  addParticle(new Particle(*newData));

  return 1; // success
}

bool ParticleTable::modifyMass(std::string name, double newMass)
{
  Particle* thisParticle = particle(name);
  if (0 != thisParticle) {
    thisParticle->data()->mass = newMass;
    return true; // success
  }
  std::cerr << "ParticleTable::modifyMass: " << name << " not found!" << std::endl;
  return false; // error
}

bool ParticleTable::modifyWidth(std::string name, double newWidth)
{
  Particle* thisParticle = particle(name);
  if (0 != thisParticle) {
    thisParticle->data()->width = newWidth;
    return true; // success
  }
  std::cerr << "ParticleTable::modifyWidth: " << name << " not found!" << std::endl;
  return false; // error
}


void ParticleTable::print(std::ostream& o) const
{
  std::map<const std::string, Particle*>::const_iterator iter = particles.begin();
  while (iter != particles.end()) {
    iter->second->print(o);
    ++iter;
  }
  return;
}

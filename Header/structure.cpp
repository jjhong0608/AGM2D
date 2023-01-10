//
// Created by 조준홍 on 2023/01/04.
//

#include "structure.h"

auto AGM::deltaFtn(double r, double h) -> double {
    return r > 2 * h ? ZEROVALUE : (UNITVALUE + std::cos(HALFVALUE * M_PI * r / h)) / (4. * h);
}

int AGM::structure::nStructures;
double AGM::structure::hf;
double AGM::structure::hs;
std::array<std::vector<AGM::structure> *, 2> *AGM::structure::structures;

AGM::structure::structure() = default;

AGM::structure::structure(int idx) : idx(idx) {}

AGM::structure::structure(int idx, const AGM::coordinate &xy) : idx(idx), xy(xy) {}

auto AGM::structure::getIdx() const -> int {
    return idx;
}

void AGM::structure::setIdx(int i) {
    structure::idx = i;
}

auto AGM::structure::getXy() const -> const AGM::coordinate & {
    return xy;
}

void AGM::structure::setXy(const AGM::coordinate &coordinate) {
    structure::xy = coordinate;
}

auto AGM::structure::getValues() const -> const AGM::value & {
    return values;
}

void AGM::structure::setValues(const AGM::value &value) {
    structure::values = value;
}

auto AGM::structure::getNStructures() -> int {
    return nStructures;
}

void AGM::structure::setNStructures(int i) {
    structure::nStructures = i;
}

auto AGM::structure::getHf() -> double {
    return hf;
}

void AGM::structure::setHf(double d) {
    structure::hf = d;
}

auto AGM::structure::getHs() -> double {
    return hs;
}

void AGM::structure::setHs(double d) {
    structure::hs = d;
}

auto AGM::structure::getStructures() -> std::array<std::vector<structure> *, 2> * {
    return structures;
}

void AGM::structure::setStructures(std::array<std::vector<structure> *, 2> *array) {
    structure::structures = array;
}

auto AGM::structure::operator-(const AGM::point &src) -> double {
    return (xy - src.getXy()).norm();
}

auto AGM::structure::operator-(const AGM::pointHeat &src) -> double {
    return (xy - src.getXy()).norm();
}

auto AGM::structure::operator-(const AGM::structure &src) -> double {
    return (xy - src.getXy()).norm();
}

auto AGM::structure::operator[](int i) -> double & {
    return xy[i];
}

auto AGM::structure::operator[](int i) const -> const double & {
    return xy[i];
}

auto AGM::structure::operator[](const std::string &string) -> double & {
    return values[string];
}

auto AGM::structure::operator[](const std::string &string) const -> const double & {
    return values[string];
}

void AGM::structure::copyFluidVelocity(AGM::pointHeat *src) {
    auto rx{std::fabs(xy[0] - src->getXy()[0])}, ry{std::fabs(xy[1] - src->getXy()[1])};
    auto dx{deltaFtn(rx, hf)}, dy{deltaFtn(ry, hf)};
    values["sol"] += src->getValue()["sol"] * dx * dy * hf * hf;
}

void AGM::structure::structureForceUpdateToFluid(AGM::pointHeat *src) {
    auto rx{std::fabs(xy[0] - src->getXy()[0])}, ry{std::fabs(xy[1] - src->getXy()[1])};
    auto dx{deltaFtn(rx, hf)}, dy{deltaFtn(ry, hf)};
    (*src)["rhs"] += values["rhs"] * dx * dy * hf * hf;
}

void AGM::structure::updateForce(std::vector<structure> *previousStructures, int component) {
    auto prev{this}, next{this};
    auto findPrevNext = [this, &prev, &next, &component]() -> void {
        if (structures->at(component)->begin()->getIdx() == getIdx()) {
            prev = &(structures->at(component)->back());
            next = &(structures->at(component)->at(1));
            return;
        } else if (structures->at(component)->back().getIdx() == getIdx()) {
            prev = &*std::prev(structures->at(component)->end() - 1);
            next = &(structures->at(component)->front());
            return;
        }
        for (auto item = structures->at(component)->begin(); item != structures->at(component)->end(); ++item) {
            if (item->getIdx() == getIdx()) {
                prev = &*std::prev(item);
                next = &*std::next(item);
                return;
            }
        }
        printError("AGM::structure::updateForce", "findPrevNext");
    };
    findPrevNext();
    auto r0_prev{previousStructures->at(getIdx()) - previousStructures->at(prev->getIdx())};
    auto r0_next{previousStructures->at(getIdx()) - previousStructures->at(next->getIdx())};
    auto ri_prev{*this - *prev};
    auto ri_next{*this - *next};
    auto d_prev{xy[component] - prev->getXy()[component]};
    auto d_next{xy[component] - next->getXy()[component]};

    values["rhs"] -= 1e2 * (ri_prev - r0_prev) * d_prev / (ri_prev * r0_prev);
    values["rhs"] -= 1e2 * (ri_next - r0_next) * d_next / (ri_next * r0_next);
}

AGM::structure::~structure() = default;

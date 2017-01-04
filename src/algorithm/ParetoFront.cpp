#include "ParetoFront.hpp"
#include "Point.hpp"
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <iostream>
#include <limits>
#include <vector>


namespace mosqp
{

ParetoFront::ParetoFront(int const max_points, size_t const num_objectives, std::vector<Point> const points)
    : maxPoints(max_points), objectiveSortings(num_objectives)
{
    AddPoints(points);
    assert(IsSortingCorrect());
}

void ParetoFront::AddPoint(Point const &new_point, bool const init)
{
    TryInsertPoint(new_point, init);
    if (points.size() > maxPoints)
    {
        Cleanup();
    }

    assert(IsSortingCorrect());
}

void ParetoFront::AddPoints(std::vector<Point> const &new_points)
{
    for (Point const &point : new_points)
    {
        TryInsertPoint(point);
    }

    if (points.size() > maxPoints)
    {
        Cleanup();
    }

    assert(IsSortingCorrect());
}

std::vector<Point>::iterator ParetoFront::RemovePoint(std::vector<Point>::const_iterator it)
{
    size_t const index = it - points.begin();
    std::vector<Point>::iterator new_it = points.erase(it);
    for (std::vector<size_t> &sorting : objectiveSortings)
    {
        for (std::vector<size_t>::iterator it = sorting.begin(); it != sorting.end(); it += 1)
        {
            if (*it == index)
            {
                sorting.erase(it);
                break;
            }
        }

        for (size_t &entry : sorting)
        {
            if (entry > index)
            {
                entry -= 1;
            }
        }

        assert(points.size() == sorting.size());
    }

    assert(IsSortingCorrect());
    return new_it;
}

bool ParetoFront::IsFull() const
{
    return NumPoints() >= maxPoints;
}

void ParetoFront::UnstopAll()
{
    for (Point &point : points)
    {
        point.SetStopped(false);
    }
}

bool ParetoFront::AllStopped() const
{
    for (Point const &point : points)
    {
        if (!point.IsStopped())
        {
            return false;
        }
    }
    return true;
}

bool ParetoFront::AllFeasible() const
{
    return (GetNumFeasible() == points.size());
}

size_t ParetoFront::GetNumFeasible() const
{
    size_t num_feasible = 0;
    for (Point const &point : points)
    {
        if (point.IsFeasible())
        {
            num_feasible += 1;
        }
    }
    return num_feasible;
}

std::vector<Point>::const_iterator ParetoFront::begin() const
{
    return points.begin();
}

std::vector<Point>::const_iterator ParetoFront::end() const
{
    return points.end();
}

std::vector<Point>::const_iterator ParetoFront::cbegin() const
{
    return points.cbegin();
}

std::vector<Point>::const_iterator ParetoFront::cend() const
{
    return points.cend();
}

void ParetoFront::TryInsertPoint(Point const &new_point, bool const init)
{
    std::vector<size_t> to_remove;
    size_t const length = points.size();

    // iterate backwards so that removing the points by index will
    // be easy afterwards
    for (size_t i = length; i > 0;)
    {
        i -= 1;
        if (!init && new_point.IsDominated(points[i]))
        {
            new_point.IsDominated(points[i]);
            return;
        }
        else if (!init && points[i].IsDominated(new_point))
        {
            to_remove.push_back(i);
        }
    }

    for (size_t i : to_remove)
    {
        RemovePoint(i);
    }

    InsertPoint(new_point);
}

void ParetoFront::InsertPoint(Point const new_point)
{
    std::vector<std::vector<size_t>::const_iterator> sorting = GetSortedIndices(new_point);
    size_t const new_point_index = points.size();
    size_t const length = objectiveSortings.size();

    points.push_back(new_point);
    for (size_t i = 0; i < length; i += 1)
    {
        objectiveSortings[i].insert(sorting[i], new_point_index);
        assert(points.size() == objectiveSortings[i].size());
    }
}

void ParetoFront::RemovePoint(size_t const index)
{
    RemovePoint(points.begin() + index);
}

void ParetoFront::Cleanup()
{
    std::vector<double> crowding_distance = ComputeCrowdingDistances();
    std::vector<double>::iterator min_element;
    while (IsOverfilled())
    {
        min_element = std::min_element(crowding_distance.begin(), crowding_distance.end());
        RemovePoint(min_element - crowding_distance.begin());
        crowding_distance.erase(min_element);
    }
}

bool ParetoFront::AllNonDominated() const
{
    size_t const num_points = NumPoints();
    for (size_t i = 0; i < num_points; i += 1)
    {
        for (size_t j = i + 1; j < num_points; j += 1)
        {
            if (points[i].IsDominated(points[j]) || points[j].IsDominated(points[i]))
            {
                return false;
            }
        }
    }

    return true;
}

std::vector<std::vector<size_t>::const_iterator> ParetoFront::GetSortedIndices(Point const &point) const
{
    size_t const num_objectives = objectiveSortings.size();
    std::vector<std::vector<size_t>::const_iterator> indices(num_objectives);
    for (size_t i = 0; i < num_objectives; i += 1)
    {
        for (indices[i] = objectiveSortings[i].begin(); indices[i] != objectiveSortings[i].end(); indices[i] += 1)
        {
            if (point.IsSmaller(points[*indices[i]], i))
            {
                break;
            }
        }
    }

    return indices;
}

std::vector<double> ParetoFront::ComputeCrowdingDistances() const
{
    std::vector<double> distances(NumPoints());
    std::fill_n(distances.begin(), NumPoints(), 0);
    size_t const num_objectives = objectiveSortings.size();
    size_t const num_points = NumPoints();
    double max_distance;

    size_t objective_index = 0;
    for (std::vector<size_t> const &objective_sorting : objectiveSortings)
    {
        distances[objective_sorting.front()] = std::numeric_limits<double>::infinity();
        distances[objective_sorting.back()] = std::numeric_limits<double>::infinity();
        max_distance = (points[objective_sorting.back()].GetObjectiveValue(objective_index) -
                        points[objective_sorting.front()].GetObjectiveValue(objective_index));
        assert(max_distance >= 0);

        for (auto it = objective_sorting.cbegin() + 1; it != objective_sorting.cend() - 1; it += 1)
        {
            distances[*it] += (points[*(it + 1)].GetObjectiveValue(objective_index) -
                               points[*(it - 1)].GetObjectiveValue(objective_index)) / max_distance;
        }

        objective_index += 1;
    }

    size_t const num_feasible = GetNumFeasible();
    if (num_feasible < maxPoints)
    {
        for (size_t i = 0; i < num_points; i += 1)
        {
            if (points[i].IsFeasible())
            {
                distances[i] = std::numeric_limits<double>::infinity();
            }
        }
    }

    return distances;
}

bool ParetoFront::IsOverfilled() const
{
    return NumPoints() > maxPoints;
}

size_t ParetoFront::NumPoints() const
{
    return points.size();
}

void ParetoFront::WriteX(std::ostream &stream)
{
    for (Point const &point : points)
    {
        for (double coordinate : point.GetX())
        {
            stream << coordinate << " ";
        }
        stream << std::endl;
    }
    stream << std::endl;
}

void ParetoFront::WriteF(std::ostream &stream)
{
    size_t const num_objectives = objectiveSortings.size();
    for (Point const &point : points)
    {
        for (size_t i = 0; i < num_objectives; i += 1)
        {
            stream << point.GetObjectiveValue(i) << " ";
        }
        stream << std::endl;
    }
    stream << std::endl;
}

bool ParetoFront::IsSortingCorrect() const
{
    size_t objective_index = 0;
    for (std::vector<size_t> const &objective_sorting : objectiveSortings)
    {
        for (int i = 1; i < objective_sorting.size(); i += 1)
        {
            if (points[objective_sorting[i]].GetObjectiveValue(objective_index) <
                points[objective_sorting[i - 1]].GetObjectiveValue(objective_index))
            {
                return false;
            }
        }
        objective_index += 1;
    }

    return true;
}

} // namespace mosqp
